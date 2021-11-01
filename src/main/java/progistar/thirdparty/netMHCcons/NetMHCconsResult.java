package progistar.thirdparty.netMHCcons;

import java.util.ArrayList;
import java.util.Hashtable;

import progistar.thirdparty.netMHCpan.HLA;

public class NetMHCconsResult {

	public String[] hlaTypes = null;
	public ArrayList<NetMHCconsData> records = new ArrayList<NetMHCconsData>();
	public Hashtable<String, NetMHCconsData> peptideToRecord = new Hashtable<String, NetMHCconsData>();
	
	public void addRecord (NetMHCconsData record) {
		this.records.add(record);
		this.peptideToRecord.put(record.peptide, record);
	}
	
	public String getHeader () {
		StringBuilder str = new StringBuilder();

		for(int i=0; i<hlaTypes.length; i++) {
			str.append(hlaTypes[i]);
			str.append("\t");
		}
		
		str.append("MHC-I\tBestType\tBestScore");
		
		return str.toString();
	}
	
	/**
	 * 
	 * Tab-delimited NB/WB/SB for each HLA type. <br>
	 * the last column is the most strong binding affinity among them. <br>
	 * 
	 * 
	 * 
	 * @param peptide
	 * @return
	 */
	public String getHLATyping (String peptide) {
		NetMHCconsData record = peptideToRecord.get(peptide);
		
		StringBuilder str = new StringBuilder();
		int maxBinding = 0;
		String bestHLAType = "-";
		String bestScore = "-";
		if(record == null) {
			for(int i=0; i<hlaTypes.length; i++) {
				str.append("NB");
				str.append("\t");
			}
		} else {
			
			bestHLAType = record.getBestHLAType();
			bestScore = record.getBestScore()+"";
			
			for(int i=0; i<hlaTypes.length; i++) {
				HLA hla = record.hlas.get(i);
				if(hla.elRank < 1250) {
					str.append("B");
					maxBinding = Math.max(1, maxBinding);
				} else {
					str.append("NB");
				}
				str.append("\t");
			}
		}
		
		if(maxBinding == 0) {
			str.append("NB\t"+bestHLAType+"\t"+bestScore);
		}else if(maxBinding == 1) {
			str.append("B\t"+bestHLAType+"\t"+bestScore);
		}
		
		return str.toString();
	}
}
