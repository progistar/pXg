package progistar.thirdparty.netMHCcons;

import java.util.ArrayList;

import progistar.thirdparty.netMHCpan.HLA;

public class NetMHCconsData {
	public String peptide;
	public ArrayList<HLA> hlas = new ArrayList<HLA>(); 
	
	public String getBestHLAType () {
		int bestHLAIndex = 0;
		
		for(int i=1; i<hlas.size(); i++) {
			if(hlas.get(i).elRank < hlas.get(bestHLAIndex).elRank) {
				bestHLAIndex = i;
			}
		}
		
		return hlas.get(bestHLAIndex).type;
	}
	
	public double getBestScore () {
		int bestHLAIndex = 0;
		
		for(int i=1; i<hlas.size(); i++) {
			if(hlas.get(i).elRank < hlas.get(bestHLAIndex).elRank) {
				bestHLAIndex = i;
			}
		}
		
		return hlas.get(bestHLAIndex).elRank;
	}
}
