package progistar.pXg.data.parser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.data.GenomicSequence;
import progistar.pXg.data.PxGAnnotation;
import progistar.pXg.data.XBlock;

public class ResultParser {

	private ResultParser () {}
	
	public static PxGAnnotation parseResult (ArrayList<File> files) {
		PxGAnnotation annotation = new PxGAnnotation();
		
		
		try {
			for(File file : files) {
				System.out.println("parsing "+file.getName()+" ...");
				BufferedReader BR = new BufferedReader(new FileReader(file));
				String line = null;
				
				String uniqueID = null;
				boolean isDecoy = false;
				double meanQScore = 0;
				
				// for unmapped reads
				String fullReads = null;
				while((line = BR.readLine()) != null) {
					String[] field = line.split("\t");
					if(field[0].equalsIgnoreCase(Constants.OUTPUT_G_UNIQUE_ID)) {
						uniqueID = field[1];
						// decoy decision
						if(uniqueID.startsWith(Constants.DECOY_PREFIX)) {
							isDecoy = true;
						} else {
							isDecoy = false;
						}
						
						fullReads = null;
						
					} else if(field[0].equalsIgnoreCase(Constants.OUTPUT_G_QSCORE)) {
						meanQScore = Double.parseDouble(field[1]);
					} else if(field[0].equalsIgnoreCase(Constants.OUTPUT_G_SEQUENCE)) {
						fullReads = field[1];
					} else if(field[0].equalsIgnoreCase(Constants.OUTPUT_G_PEPTIDE)) {
						
						try {
							String pSeq = field[1]; // peptide sequence without I/L consideration
							XBlock xBlock = new XBlock();
							xBlock.genomicLocus = field[2];
							xBlock.strand = field[3].charAt(0);
							
							xBlock.leftFlankSequence = field[4];
							xBlock.genomicSequence = field[5];
							xBlock.rightFlankSequence = field[6];
							
							xBlock.leftFlankRefSequence = field[7];
							xBlock.referenceSequence = field[8];
							xBlock.rightFlankRefSequence = field[9];
							
							xBlock.mutations = field[10];
							xBlock.mutationStatus = field[11];
							xBlock.tAnnotations = field[12];
							
							String[] dists = field[13].split("\\|");
							xBlock.fullReadSequence = fullReads;
							// If unmapped reads, merging xBlocks and making a single contig xBlock.
							
							// antisense checker if three-frame
							/**
							 * Antisense event is not allowed.
							 * TODO: More events can be added in future.
							 */
							if(Parameters.translationMethod == Constants.THREE_FRAME) {
								StringBuilder transcriptsWithOutBANlist = new StringBuilder();
								StringBuilder distsWithOutBANlist = new StringBuilder();
								
								String[] transcripts = xBlock.tAnnotations.split("\\|");
								for(int i=0; i<transcripts.length; i++) {
									String transcript = transcripts[i];
									if(!transcript.contains(";antisense;")) {
										if(transcriptsWithOutBANlist.length() != 0) {
											transcriptsWithOutBANlist.append("|");
											distsWithOutBANlist.append("|");
										}
										transcriptsWithOutBANlist.append(transcript);
										distsWithOutBANlist.append(dists[i]);
									}
								}
								
								// If all transcripts were discarded, then it means only possible interpretation is antisense
								// in this case, we accept antisense even three-frame translation.
								if(transcriptsWithOutBANlist.length() != 0) {
									xBlock.tAnnotations = transcriptsWithOutBANlist.toString();
									field[13] = distsWithOutBANlist.toString();
								}
							}
							
							// calculate dists (longest exon length)
							dists = field[13].split("\\|");
							String percentFullDist = "-";
							String percentExonDist = "-";
							
							int longestValue = 0;
							for(int i=0; i<dists.length; i++) {
								String[] values = dists[i].split("\\;");
								if(!values[1].equalsIgnoreCase("-")) {
									String exonLength = values[1].split("\\/")[1];
									int value = Integer.parseInt(exonLength);
									
									if(value > longestValue) {
										longestValue = value;
										percentFullDist = values[0];
										percentExonDist = values[1];
									}
								}
							}
							
							if(!percentFullDist.equalsIgnoreCase("-")) {
								double v1 = Double.parseDouble(percentFullDist.split("\\/")[0]);
								double v2 = Double.parseDouble(percentFullDist.split("\\/")[1]);
								percentFullDist = (v1/v2)+"";
							}
							
							if(!percentExonDist.equalsIgnoreCase("-")) {
								double v1 = Double.parseDouble(percentExonDist.split("\\/")[0]);
								double v2 = Double.parseDouble(percentExonDist.split("\\/")[1]);
								percentExonDist = (v1/v2)+"";
							}
							
							xBlock.percentFullDist = percentFullDist;
							xBlock.percentExonDist = percentExonDist;
							
							if(xBlock.strand == '+') {
								xBlock.peptideSequence = GenomicSequence.translation(xBlock.genomicSequence, 0);
							} else {
								xBlock.peptideSequence = GenomicSequence.reverseComplementTranslation(xBlock.genomicSequence, 0);
							}
							
							if(xBlock.mockReadCount != 0 || xBlock.targetReadCount != 0 || xBlock.decoyQScore != 0 || xBlock.targetQScore != 0) {
								System.out.println("Line 100 at ResultParser.java: why they are not zero?");
							}
							
							if(isDecoy) {
								xBlock.mockReadCount++;
								xBlock.decoyQScore += meanQScore;
							} else {
								xBlock.targetReadCount++;
								xBlock.targetQScore += meanQScore;
							}
							
							xBlock.sequenceID = uniqueID;
							if(xBlock.tAnnotations.contains(Constants.MARK_UNMAPPED+"")) {
								// those xBlocks are needed to merge into contig.
								xBlock.genomicLocus = "-";
								// TODO: assembly?
								// Assembler.addXBlock(pSeq, xBlock);
							}
							
							annotation.putXBlock(pSeq, xBlock);
						}catch(Exception e) {
							e.printStackTrace();
							System.out.println(line);
						}
					}
				}
				
				BR.close();
			}
		}catch(IOException ioe ) {
			
		}
		
		// just write down unmapped reads (which mapped to any peptides)
		//Assembler.write();
		
		/*
		// assemble xBlocks
		ArrayList<XBlock> xBlocks = Assembler.assemble();
		
		xBlocks.forEach(xBlock -> {
			String pSeq = Assembler.getPSeq(xBlock);
			annotation.putXBlock(pSeq, xBlock);
		});
		*/
		

		// TODO: v2.0.0
		// NetMHCpan option!
		if(true) {
			
		}
		
		
		return annotation;
	}
}
