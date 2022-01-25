package progistar.thirdparty.netMHCpan;

import java.util.ArrayList;
import java.util.Hashtable;

public class NetMHCpanResult {

	public String[] hlaTypes = null;
	public ArrayList<NetMHCpanData> records = new ArrayList<NetMHCpanData>();
	public Hashtable<String, NetMHCpanData> peptideToRecord = new Hashtable<String, NetMHCpanData>();
	
	public void addRecord (NetMHCpanData record) {
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
		NetMHCpanData record = peptideToRecord.get(peptide);
		
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
				
				if(hla.elRank < 0.5) {
					str.append(hla.elRank+"");
					maxBinding = Math.max(2, maxBinding);
				} else if(hla.elRank < 2) {
					str.append(hla.elRank+"");
					maxBinding = Math.max(1, maxBinding);
				} else {
					str.append(hla.elRank+"");
				}
				str.append("\t");
			}
		}
		
		if(maxBinding == 0) {
			str.append("NB").append("\t"+bestHLAType+"\t"+bestScore);
		}else if(maxBinding == 1) {
			str.append("WB").append("\t"+bestHLAType+"\t"+bestScore);
		}else if(maxBinding == 2) {
			str.append("SB").append("\t"+bestHLAType+"\t"+bestScore);
		}
		
		return str.toString();
	}
}
