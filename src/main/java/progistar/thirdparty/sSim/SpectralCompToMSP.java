package progistar.thirdparty.sSim;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import progistar.thirdparty.netMHCpan.NetMHCpanParser;
import progistar.thirdparty.netMHCpan.NetMHCpanResult;

public class SpectralCompToMSP {

	public static void main(String[] args) throws IOException {
		String[] exSpectraSet = {
				//"/Users/gistar/projects/pXg/MSMS/B_LCL1.mgf",
				//"/Users/gistar/projects/pXg/MSMS/B_LCL2.mgf",
				//"/Users/gistar/projects/pXg/MSMS/B_LCL3.mgf",
				//"/Users/gistar/projects/pXg/MSMS/B_LCL4.mgf",
				"/Users/gistar/projects/pXg/MSMS/DOHH2_400M_050219.mgf",
				"/Users/gistar/projects/pXg/MSMS/HBL1_DMSO_200M_050219.mgf",
				"/Users/gistar/projects/pXg/MSMS/SUDHL4_400M_050219.mgf"
				//"/Users/gistar/projects/pXg/MSMS/THP1_1.mgf",
				//"/Users/gistar/projects/pXg/MSMS/THP1_2.mgf",
				//"/Users/gistar/projects/pXg/MSMS/THP1_3.mgf"
				
		};
		
		String[] predSpectraSet = {
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_fullreads/DOHH2.msp",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_fullreads/HBL1.msp",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_fullreads/SUDHL4.msp"
		};
		
		String[] deepLCSet = {
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_fullreads/PEAKS_DOHH2.nocut.input_deeplc_predictions.csv",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_fullreads/PEAKS_HBL1.nocut.input_deeplc_predictions.csv",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_fullreads/PEAKS_SUDHL4.nocut.input_deeplc_predictions.csv"
		};
		
		String[] netmhcpanSet = {
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_fullreads/DOHH2.netmhcpan.xls",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_fullreads/HBL1.netmhcpan.xls",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_fullreads/SUDHL4.netmhcpan.xls"
		};
		
		String[] pinSet = {
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_fullreads/PEAKS_DOHH2.nocut.pxg.pin",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_fullreads/PEAKS_HBL1.nocut.pxg.pin",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_fullreads/PEAKS_SUDHL4.nocut.pxg.pin"
				
		};
		
		
		for(int idx = 0; idx<pinSet.length; idx++) {
			System.out.println(pinSet[idx]);
			// index by title
			Spectra exSpectra = new Spectra(exSpectraSet[idx], Spectra.FILE_TYPE_MGF);
			// index by peptide/charge
			Spectra predSpectra = new Spectra(predSpectraSet[idx], Spectra.FILE_TYPE_MSP);
			// deepLC best
			File deepLCRes = new File(deepLCSet[idx]);
			// HLA binding prediction
			File netmhcpan = new File(netmhcpanSet[idx]);
			
			NetMHCpanResult netMHCpanResult = NetMHCpanParser.parseNetMHCpan(netmhcpan.getAbsolutePath());
			
			File pinFile = new File(pinSet[idx]);
			String line = null;
			
			BufferedReader BRlc = new BufferedReader(new FileReader(deepLCRes));
			Hashtable<String, Double> bestDelta = new Hashtable<String, Double>();
			Hashtable<String, ArrayList<Double>> deltaRTs = new Hashtable<String, ArrayList<Double>>();
			BRlc.readLine();// skip header
			
			while((line = BRlc.readLine()) != null) {
				String[] fields = line.split("\\,");
				String peptide =fields[1];
				Double obRT = Double.parseDouble(fields[3]);
				Double predRT = Double.parseDouble(fields[4]);
				double delta = obRT - predRT;
				Double thisDelta = bestDelta.get(peptide);
				ArrayList<Double> deltaList = deltaRTs.get(peptide);
				if(thisDelta == null || Math.abs(thisDelta) > Math.abs(delta)) {
					bestDelta.put(peptide, delta);
				}
				
				if(deltaList == null) {
					deltaList = new ArrayList<Double>();
				}
				deltaList.add(delta);
				deltaRTs.put(peptide, deltaList);
			}
			
			BRlc.close();
			
			BufferedReader BR = new BufferedReader(new FileReader(pinFile));
			BufferedWriter BW = new BufferedWriter(new FileWriter(pinFile.getAbsolutePath().replace(".pin", ".predfeat.pin")));
			
			
			StringBuilder output = new StringBuilder();
			String[] headers = BR.readLine().split("\t");
			for(int i=0; i<headers.length-2; i++) {
				output.append(headers[i]).append("\t");
			}
			output.append("SA\tBestDeltaRT\tAvgDelta\tmLog2BestELRank\t").append(headers[headers.length-2]).append("\t").append(headers[headers.length-1]);
			BW.append(output.toString());
			BW.newLine();
			while((line = BR.readLine()) != null) {
				output.setLength(0);
				String[] fields = line.split("\t");
				String uniqueID = fields[0];
				String fileName = uniqueID.split("\\|")[0].replace(".raw", "");
				String scanNum = uniqueID.split("\\|")[1].split("\\:")[1];
				String charge = uniqueID.split("\\|")[2];
				String title = fileName+"."+scanNum+"."+scanNum+"."+charge;
				String peptide = fields[fields.length-2];
				
				Spectrum s1 = exSpectra.getSpectrumByScanNum(title);
				Spectrum s2 = predSpectra.getSpectrumByScanNum(peptide+"/"+charge);
				
				s1.setPeptide(peptide);
				s2.setPeptide(peptide);
				
				double scaScore = SpectralScores.spectralContrastAngle(s1, s2, 0.02, false);
				
				for(int i=0; i<fields.length-2; i++) {
					output.append(fields[i]).append("\t");
				}
				
				ArrayList<Double> deltaList = deltaRTs.get(peptide);
				double avgDelta = 0;
				for(Double delta : deltaList) {
					avgDelta += delta;
				}
				avgDelta /= (double) deltaList.size();
				
				double elrank = -Math.log(netMHCpanResult.peptideToRecord.get(peptide).getBestScore()) / Math.log(2);
				
				output.append(scaScore).append("\t")
				.append(bestDelta.get(peptide)).append("\t")
				.append(avgDelta).append("\t")
				.append(elrank).append("\t")
				.append(fields[fields.length-2]).append("\t").append(fields[fields.length-1]);
				BW.append(output.toString());
				BW.newLine();
				
			}
			BW.close();
			BR.close();
		}
		
	}
	
	
}
