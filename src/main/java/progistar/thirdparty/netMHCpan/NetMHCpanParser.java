package progistar.thirdparty.netMHCpan;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class NetMHCpanParser {

	public static NetMHCpanResult parseNetMHCpan (String fileName) {
		
		NetMHCpanResult result = new NetMHCpanResult();
		
		try {
			BufferedReader BR = new BufferedReader(new FileReader(fileName));
			
			String[] headers = BR.readLine().split("\\s");
			int sizeOfHLAs = 0;
			for(int i=0; i<headers.length; i++) {
				if(headers[i].length() != 0) {
					sizeOfHLAs++;
				}
			}
			String[] hlaTypes = new String[sizeOfHLAs];
			sizeOfHLAs = 0;
			for(int i=0; i<headers.length; i++) {
				if(headers[i].length() != 0) {
					hlaTypes[sizeOfHLAs++] = headers[i];
				}
			}
			
			result.hlaTypes = hlaTypes;

			int peptideIndex = 1;
			int elRankBaseIndex = 6;
			int elRankIndexInteravl = 4;
			
			String line = BR.readLine(); // skip header
			
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String peptide = fields[peptideIndex];
				NetMHCpanData data = new NetMHCpanData();
				
				data.peptide = peptide;
				
				double[] elRanks = new double[hlaTypes.length];
				
				for(int i=0; i<elRanks.length; i++) {
					elRanks[i] = Double.parseDouble(fields[elRankBaseIndex+elRankIndexInteravl*i]);
					
					HLA hla = new HLA();
					hla.elRank = elRanks[i];
					hla.type = hlaTypes[i];
					
					data.hlas.add(hla);
				}
				
				result.addRecord(data);
			}
			
			BR.close();
		}catch(IOException ioe) {
			
		}
		
		System.out.println("a total of "+result.records.size()+" was inserted..");
		
		return result;
	}
}
