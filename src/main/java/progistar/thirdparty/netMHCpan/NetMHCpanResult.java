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
		
		str.append("MHC-I").append("\t").append("elRank");
		
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
		double bestBindingAffinity = 1;
		if(record == null) {
			for(int i=0; i<hlaTypes.length; i++) {
				str.append("NB");
				str.append("\t");
			}
		} else {
			for(int i=0; i<hlaTypes.length; i++) {
				HLA hla = record.hlas.get(i);
				
				bestBindingAffinity = Math.min(bestBindingAffinity, hla.elRank/100);
				
				if(hla.elRank < 0.5) {
					str.append("SB");
					maxBinding = Math.max(1, maxBinding);
				} else if(hla.elRank < 2) {
					str.append("WB");
					maxBinding = Math.max(2, maxBinding);
				} else {
					str.append("NB");
				}
				str.append("\t");
			}
		}
		
		if(maxBinding == 0) {
			str.append("NB").append("\t"+bestBindingAffinity);
		}else if(maxBinding == 1) {
			str.append("WB").append("\t"+bestBindingAffinity);
		}else if(maxBinding == 2) {
			str.append("SB").append("\t"+bestBindingAffinity);
		}
		
		return str.toString();
	}
}
