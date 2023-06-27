package progistar.pXg.data.parser;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

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
		BufferedReader BR = new BufferedReader(new FileReader(Parameters.sequenceFilePath));
		String line = null;
		
		while((line = BR.readLine()) != null) {
			if(line.startsWith("@")) {
				BW.append(line);
				BW.newLine();
			} else {
				String sequenceID = line.split("\\s")[0];
				// for target check
				if(sequenceIDChecker.get(sequenceID) != null) {
					BW.append(line);
					BW.newLine();
				}
			}
		}
		
		BR.close();
	}
}
