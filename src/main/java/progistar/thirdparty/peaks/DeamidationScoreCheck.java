package progistar.thirdparty.peaks;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

public class DeamidationScoreCheck {

	static class PEAKSRes {
		public static final int FILE_IDX = 1;
		public static final int PEPT_IDX = 3;
		public static final int SCAN_IDX = 4;
		public static final int ALC_IDX = 7;
		public static final int ALC_SINGLE_IDX = 15;
		
		public int deltaScore = 100;
		public boolean isDeamidatedAssociated = false;
		public String counterpart = null;
		
		// deamiNumber of top-ranked > deamiNumber of other-ranked -> 1
		// deamiNumber of top-ranked = deamiNumber of other-ranked -> 0
		// deamiNumber of top-ranked < deamiNumber of other-ranked -> -1
		public int category = 0;
		
		public ArrayList<String> candidates = new ArrayList<String>();
		
		public void calALC () {
			for(int i=0; i<candidates.size(); i++) {
				double alc = 0;
				String[] fields = candidates.get(i).split("\t");
				
				String[] localLCs = fields[ALC_SINGLE_IDX].split("\\s");
				for(String lc : localLCs) {
					alc += Double.parseDouble(lc);
				}
				
				alc /= localLCs.length;
				
				int alcVal = Integer.parseInt(fields[ALC_IDX]);
				if(alcVal != Math.floor(alc)) {
					System.out.println(alcVal+"\t"+alc+"\t"+fields[ALC_SINGLE_IDX]);
				}
			}
			
		}
		
		public void compareDeamiIDs () {
			String[] fields = candidates.get(0).split("\t");
			String peptide = fields[PEPT_IDX];
			int score = Integer.parseInt(fields[ALC_IDX]);
			if(peptide.contains("N") || peptide.contains("Q") || peptide.contains("D") || peptide.contains("E")) {
				String altPeptide = peptide.replace("N(+.98)", "D");
				altPeptide = altPeptide.replace(("Q(+.98)"), "E");
				
				for(int i=1; i<candidates.size(); i++) {
					String[] thisFields = candidates.get(i).split("\t");
					String thisPeptide = thisFields[PEPT_IDX];
					thisPeptide = thisPeptide.replace("N(+.98)", "D");
					thisPeptide = thisPeptide.replace("Q(+.98)", "E");
					
					if(altPeptide.equalsIgnoreCase(thisPeptide)) {
						this.deltaScore = score - Integer.parseInt(thisFields[ALC_IDX]);
						counterpart = candidates.get(i);
						isDeamidatedAssociated = true;
						
						int deamiCount = 0;
						int thisDeamiCount = 0;
						
						altPeptide = peptide;
						thisPeptide = thisFields[PEPT_IDX];
						
						while(altPeptide.contains(".98")) {
							altPeptide = altPeptide.replaceFirst(".98", "");
							deamiCount ++;
						}
						
						while(thisPeptide.contains(".98")) {
							thisPeptide = thisPeptide.replaceFirst(".98", "");
							thisDeamiCount ++;
						}
						
						if(deamiCount > thisDeamiCount) {
							category = 1;
						} else if(deamiCount < thisDeamiCount) {
							category = -1;
						} else {
							category = 0;
						}
						
						break;
					}
				}
			}
		}
		
		
	}
	
	public static void main(String[] args) throws IOException {
		File peaksAllResFile = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/3.withCalibrationAddScanNum/pXg/PeptideAnnotationS1_5ppm_002_MSFrecal_addScanNum.tsv");
		
		BufferedReader BR = new BufferedReader(new FileReader(peaksAllResFile));
		String line = null;
		
		BR.readLine(); // skip header
		
		Hashtable<String, PEAKSRes> peaksRes = new Hashtable<String, PEAKSRes>();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String scanID = fields[PEAKSRes.FILE_IDX] +"_"+fields[PEAKSRes.SCAN_IDX];
			
			PEAKSRes record = peaksRes.get(scanID);
			if(record == null) {
				record = new PEAKSRes();
			}
			
			record.candidates.add(line);
			peaksRes.put(scanID,record);
		}
		
		BR.close();
		
		System.out.println("Let's Compare");
		System.out.println("DeltaALCScore,ScanID,Top-ranked,Counterpart,Class");
		peaksRes.forEach((scanID, res) -> {
			res.compareDeamiIDs();
			
			if(res.isDeamidatedAssociated) {
				int category = res.category;
				String categoryStr = "same";
				if(category == 1) {
					categoryStr = "greater";
				} else if(category == -1) {
					categoryStr = "less";
				}
				res.calALC();
//				System.out.println(res.deltaScore+","+scanID+","+res.candidates.get(0).split("\t")[PEAKSRes.PEPT_IDX]+","+res.counterpart.split("\t")[PEAKSRes.PEPT_IDX]+","+categoryStr);
			}
			
//					if(res.deltaScores[i] >= 2) {
//						System.out.println(res.candidates.get(i));
//						System.out.println(res.compLines[i]);
//						System.out.println("------");
//					}
		});
		
		
	}
}
