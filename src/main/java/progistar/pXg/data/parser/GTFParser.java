package progistar.pXg.data.parser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import progistar.pXg.constants.Constants;
import progistar.pXg.data.ABlock;
import progistar.pXg.data.GenomicAnnotation;
import progistar.pXg.data.TBlock;
import progistar.pXg.utils.IndexConvertor;

public class GTFParser {

	// prevent to generate constructor
	private GTFParser () {}

	private static final String[] SELECTED_FEATURES = {"exon", "cds"};


	private enum FieldIndex {
		CHR(0), SOURCE(1), FEATURE(2), START(3), END(4), STRAND(6), ATTR(8);

		private int value;

		FieldIndex(int value) {
			this.value = value;
		}
	}

	private static String getGtfAttr(String[] attr, String tag){

		for(String _s : attr){
			if(_s.contains(tag)){
				return _s.replaceAll("[\"\\s]|"+tag, "");
			}
		}

		return null;
	}

	public static GenomicAnnotation parseGTF (String gtfFilePath) {
		System.out.print("Parsing GTF: "+gtfFilePath);
		long startTime = System.currentTimeMillis();

		GenomicAnnotation annotation = new GenomicAnnotation();
		try {
			File samFile = new File(gtfFilePath);

			BufferedReader BR = new BufferedReader(new FileReader(samFile));

			String line = null;

			int chrIndex = FieldIndex.CHR.value;
			int startIndex = FieldIndex.START.value;
			int endIndex = FieldIndex.END.value;
			int featureIndex = FieldIndex.FEATURE.value;
			int strandIndex = FieldIndex.STRAND.value;
			int attrIndex = FieldIndex.ATTR.value;

			// add transcript information and corresponding exon structures
			while((line = BR.readLine()) != null) {
				if(line.startsWith("#"))
				 {
					continue; // skip meta
				}

				String[] fields = line.split("\t");

				String feature = fields[featureIndex];
				int start = Integer.parseInt(fields[startIndex]);
				int end = Integer.parseInt(fields[endIndex]);
				String[] attr = fields[attrIndex].split(";");

				if(feature.equalsIgnoreCase("transcript")) {
					String transcriptID = getGtfAttr(attr, "transcript_id");
					String transcriptName = getGtfAttr(attr, "transcript_name");
					String transcriptType = getGtfAttr(attr, "transcript_type");

					String geneID = getGtfAttr(attr, "gene_id");
					String geneName = getGtfAttr(attr, "gene_name");
					String geneType = getGtfAttr(attr, "gene_type");

					boolean strand = fields[strandIndex].equalsIgnoreCase("-") ? false : true;

					String chr = fields[chrIndex];
					// enroll chr index and chr string
					IndexConvertor.putChrIndexer(chr);

					// put transcript ID
					annotation.putTBlock(IndexConvertor.chrToIndex(chr), strand, start, end,
							transcriptID, transcriptName, transcriptType, geneID, geneName, geneType);
				}

				else {
					boolean isSelectedFeature = false;
					for(String sFeature : SELECTED_FEATURES) {
						if(feature.equalsIgnoreCase(sFeature)) {
							isSelectedFeature = true;
						}
					}

					// skip if it is not a selected feature.
					if(!isSelectedFeature) {
						continue;
					}

					// selected features consist of structural blocks which are building block of a transcript.
					byte bFeature = feature.equalsIgnoreCase("CDS") ? Constants.CDS : Constants.EXON;


					String transcriptID = getGtfAttr(attr, "transcript_id");
					TBlock tBlock = annotation.getTBlockByTID(transcriptID);

					// ASSERT!
					assert tBlock != null;

					ABlock aBlock = new ABlock();
					aBlock.start = start;
					aBlock.end = end;
					aBlock.feature = bFeature;
					aBlock.transcriptIndex = tBlock.tBlockID;

					tBlock.aBlocks.add(aBlock);
				}
			}

			BR.close();

			// classify exon into CDS, NCDS and UTR.
			// plus, define intron and intergenic regions.
			annotation.assignTypesInTBlocks();

		}catch(IOException ioe) {
			System.out.println("...\tFail to load GTF");
		}


		// update ENST mapper
		annotation.updateENSTMapper();

		long endTime = System.currentTimeMillis();
		System.out.println("\tElapsed time: "+((endTime-startTime)/1000) + " sec");

		return annotation;
	}
}
