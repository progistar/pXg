package progistar.pXg.data.parser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import progistar.pXg.constants.Constants;
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
				
				while((line = BR.readLine()) != null) {
					String[] field = line.split("\t");
					if(field[0].equalsIgnoreCase(Constants.OUTPUT_G_UNIQUE_ID)) {
						uniqueID = field[1];
						// decoy decision
						if(uniqueID.startsWith("XXX")) {
							isDecoy = true;
						} else {
							isDecoy = false;
						}
					} else if(field[0].equalsIgnoreCase(Constants.OUTPUT_G_PEPTIDE)) {
						
						String pSeq = field[1]; // peptide sequence without I/L consideration
						XBlock xBlock = new XBlock();
						xBlock.genomicLocus = field[2];
						xBlock.strand = field[3].charAt(0);
						xBlock.genomicSequence = field[4];
						xBlock.tAnnotations = field[5];
						if(xBlock.strand == '+') {
							xBlock.peptideSequence = GenomicSequence.translation(xBlock.genomicSequence, 0);
						} else {
							xBlock.peptideSequence = GenomicSequence.reverseComplementTranslation(xBlock.genomicSequence, 0);
						}
						
						if(isDecoy) {
							xBlock.decoyReadCount++;
						} else {
							xBlock.targetReadCount++;
						}
						
						// put xBlock
						annotation.putXBlock(pSeq, xBlock);
					}
				}
				
				BR.close();
			}
		}catch(IOException ioe ) {
			
		}
		
		
		return annotation;
	}
}
