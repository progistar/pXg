package progistar.pXg.data.parser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Hashtable;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.data.XBlock;

public class SAMExportor {

	private static Hashtable<String, String> sequenceIDChecker = new Hashtable<String, String>();
	
	public static void putSequenceID (XBlock xBlock) {
		
		sequenceIDChecker.put(xBlock.sequenceID.split("\\@")[0].replace(Constants.DECOY_PREFIX, ""), "");
		
		xBlock.siblingXBlocks.forEach((sxBlock) -> {
			sequenceIDChecker.put(sxBlock.sequenceID.split("\\@")[0].replace(Constants.DECOY_PREFIX, ""), "");
		});
		
	}
	
	public static void writeSAM (BufferedWriter BW) throws IOException {
		
		try (SamReader samReader = SamReaderFactory.makeDefault().open(new File(Parameters.sequenceFilePaths[Parameters.CURRENT_FILE_INDEX]))) {
			
			String lineSeparator = System.getProperty("line.separator");
			
        	// get records
            for (SAMRecord samRecord : samReader) {
                // Process each SAMRecord as needed
            	String record = samRecord.getSAMString().replace(lineSeparator, "");
            	String sequenceID = record.split("\\s")[0];
				// for target check
				if(sequenceIDChecker.get(sequenceID) != null) {
					BW.append(record);
					BW.newLine();
				}
            }
            
            
        } catch (IOException e) {
            e.printStackTrace();
        }
	}
}
