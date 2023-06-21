package progistar.thirdparty.sSim;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.Iterator;

class SCA {
	String peptide;
	Spectrum spectrum1;
	Spectrum spectrum2;
	double score;
}

public class SpectralScores {

	// KVGAVVHLK_3
	// test
	public static void main(String[] args) throws IOException {
		
		
		
		
		PeptideLoader.loadPeptideList("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/ProteomeToolsHLA/Matched_NCPeptideList.tsv");
		
		// load pxg MGF
		Spectra pxgSpectra = new Spectra("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/ProteomeToolsHLA/selectedMGF_pXg.mgf", Spectra.FILE_TYPE_MGF);
		Spectra ptSpectra = new Spectra("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/ProteomeToolsHLA/selectedMGF_PT.mgf", Spectra.FILE_TYPE_MGF);
		
		calSCA(pxgSpectra, ptSpectra, pxgSpectra);
	}
	
	public static ArrayList<SCA> calSCA (Spectra spectra1, Spectra spectra2, Spectra filterSpectra) throws IOException {
		ArrayList<SCA> scas = new ArrayList<SCA>();
		
		// aggregate spectra1 by charge and peptide
		Hashtable<String, ArrayList<Spectrum>> aggregatedSpec1 = aggregateSpectra(spectra1);
		Hashtable<String, ArrayList<Spectrum>> aggregatedSpec2 = aggregateSpectra(spectra2);
		Hashtable<String, ArrayList<Spectrum>> aggregatedSpecFilter = aggregateSpectra(filterSpectra);
		
		Iterator<String> keys = (Iterator<String>) aggregatedSpecFilter.keys();
		BufferedWriter BW = new BufferedWriter(new FileWriter("matchedMGF.mgf"));
		while(keys.hasNext()) {
			String key = keys.next();
			ArrayList<Spectrum> s1 = aggregatedSpec1.get(key);
			ArrayList<Spectrum> s2 = aggregatedSpec2.get(key);
			
			if(s1 == null || s2 == null) {
				
			} else {
				Spectrum maxS1 = null;
				Spectrum maxS2 = null;
				double maxSCAScore = 0;
				
				for(Spectrum s1_ : s1) {
					for(Spectrum s2_ : s2) {
						
						if(s1_.getTitle().equalsIgnoreCase(s2_.getTitle())) {
							continue;
						}
						
						
						double scaScore = spectralContrastAngle(s1_, s2_, 0.02, false);
						
						if(scaScore > maxSCAScore) {
							maxS1 = s1_;
							maxS2 = s2_;
							maxSCAScore = scaScore;
						}
//						System.out.println(scaScore);
					}
				}
				
				if(maxS1 != null && maxS2 != null ) {
					String[] targetFileNameSplit = maxS1.getTitle().split("\\/");
					String[] pivotFileNameSplit = maxS2.getTitle().split("\\/");
					System.out.println(key+"\t"+s1.size()+"\t"+s2.size()+"\t"+targetFileNameSplit[targetFileNameSplit.length-1]+
							"\t"+pivotFileNameSplit[pivotFileNameSplit.length-1]+"\t"+maxSCAScore);
					
					BW.append(key+"\t"+maxS1.getTitle());
					BW.newLine();
					for(double[] peak : maxS1.peaks) {
						BW.append(peak[0]+"\t"+String.format("%.12f", peak[1]));
						BW.newLine();
					}
					
					BW.append(key+"\t"+maxS2.getTitle());
					BW.newLine();
					for(double[] peak : maxS2.peaks) {
						BW.append(peak[0]+"\t"+String.format("%.12f", peak[1]));
						BW.newLine();
					}
					
				}
				
			}
		}
		BW.close();
		
		
		return scas;
	}
	
	public static Hashtable<String, ArrayList<Spectrum>> aggregateSpectra (Spectra spectra) {
		Hashtable<String, ArrayList<Spectrum>> spectraPerPrecursor = new Hashtable<String, ArrayList<Spectrum>>();
		
		for(int i=0; i<spectra.sizeOfSpectra(); i++) {
			Spectrum s = spectra.getSpectrumByIndex(i);
			String peptide = s.getPeptide();
			int charge = s.getCharge();
			
			String key = peptide+"_"+charge;
			ArrayList<Spectrum> selectedSpectra = spectraPerPrecursor.get(key);
			if(selectedSpectra == null) {
				selectedSpectra = new ArrayList<Spectrum>();
				spectraPerPrecursor.put(key, selectedSpectra);
			}
			selectedSpectra.add(s);
		}
		
		System.out.println("Aggregate spectra... :" + spectra.sizeOfSpectra()+"=>"+spectraPerPrecursor.size());
		return spectraPerPrecursor;
				
	}
	
	public static ArrayList<double[]> getAnnotatedPeaks(Spectrum s1_, double tolerance) {
		String peptide = s1_.getPeptide();
		Peptide p = new Peptide(peptide, "", "");
		ArrayList<double[]> precursorIonS1 = new ArrayList<double[]>();
		precursorIonS1.addAll(s1_.getSubPeaks(s1_.getPrecursorMz()-tolerance, s1_.getPrecursorMz()+tolerance));
		precursorIonS1.addAll(s1_.getSubPeaks(s1_.getPrecursorMz()-(ProteomeConstants.H2O/((double)s1_.getCharge()))-tolerance, s1_.getPrecursorMz()-(ProteomeConstants.H2O/((double)s1_.getCharge()))+tolerance));
		precursorIonS1.addAll(s1_.getSubPeaks(s1_.getPrecursorMz()-(ProteomeConstants.NH3/((double)s1_.getCharge()))-tolerance, s1_.getPrecursorMz()-(ProteomeConstants.NH3/((double)s1_.getCharge()))+tolerance));
		
		ArrayList<double[]> experimentalPeaks = new ArrayList<double[]>();
		
		
		for(int c=1; c<s1_.getCharge(); c++) {
			double[] peaks = p.getTheoreticalLadder(ProteomeConstants.Y_ION, c, false);
			for(int i=0; i<peaks.length; i++) {
				experimentalPeaks.addAll(s1_.getSubPeaks(peaks[i]-0.02, peaks[i]+0.02));
			}
			
			peaks = p.getTheoreticalLadder(ProteomeConstants.B_ION, c, false);
			for(int i=0; i<peaks.length; i++) {
				experimentalPeaks.addAll(s1_.getSubPeaks(peaks[i]-0.02, peaks[i]+0.02));
			}
		}
		
		if(experimentalPeaks.size() != 0) {
			Collections.sort(experimentalPeaks, new Comparator<double[]>() {
				@Override
				public int compare(double[] o1, double[] o2) {
					if(o1[0] < o2[0]) {
						return -1;
					}else if(o1[0] > o2[0]) {
						return 1;
					}
					return 0;
				}
				
			});
		}
		return experimentalPeaks;
	}
	
	/**
	 * experimentalPeaks and theoreticalPeaks must be sorted by ascending order of mz. <br> 
	 * Tolerance is da. <br>
	 * 
	 * @param onePeaks
	 * @param otherPeaks
	 * @param tolerance
	 * @return
	 */
	public static double spectralContrastAngle (Spectrum s1_, Spectrum s2_, double tolerance, boolean isSquareRootNorm) {
		
		// max SI
		double maxExp = 1;
		double maxThr = 1;
		String peptide = s1_.getPeptide();
		Peptide p = new Peptide(peptide, "", "");
		
		ArrayList<double[]> precursorIonS1 = new ArrayList<double[]>();
		ArrayList<double[]> precursorIonS2 = new ArrayList<double[]>();
		
		precursorIonS1.addAll(s1_.getSubPeaks(s1_.getPrecursorMz()-tolerance, s1_.getPrecursorMz()+tolerance));
		precursorIonS1.addAll(s1_.getSubPeaks(s1_.getPrecursorMz()-(ProteomeConstants.H2O/((double)s1_.getCharge()))-tolerance, s1_.getPrecursorMz()-(ProteomeConstants.H2O/((double)s1_.getCharge()))+tolerance));
		precursorIonS1.addAll(s1_.getSubPeaks(s1_.getPrecursorMz()-(ProteomeConstants.NH3/((double)s1_.getCharge()))-tolerance, s1_.getPrecursorMz()-(ProteomeConstants.NH3/((double)s1_.getCharge()))+tolerance));
		
		precursorIonS2.addAll(s2_.getSubPeaks(s2_.getPrecursorMz()-tolerance, s2_.getPrecursorMz()+tolerance));
		precursorIonS2.addAll(s2_.getSubPeaks(s2_.getPrecursorMz()-(ProteomeConstants.H2O/((double)s2_.getCharge()))-tolerance, s2_.getPrecursorMz()-(ProteomeConstants.H2O/((double)s2_.getCharge()))+tolerance));
		precursorIonS2.addAll(s2_.getSubPeaks(s2_.getPrecursorMz()-(ProteomeConstants.NH3/((double)s2_.getCharge()))-tolerance, s2_.getPrecursorMz()-(ProteomeConstants.NH3/((double)s2_.getCharge()))+tolerance));
		
		for(double[] peak : precursorIonS1) {
			for(int i=0; i<s1_.sizeOfPeaks(); i++) {
				if(peak[0] == s1_.getPeak(i)[0]) {
					s1_.removePeak(i--);
				}
			}
		}
		
		for(double[] peak : precursorIonS2) {
			for(int i=0; i<s2_.sizeOfPeaks(); i++) {
				if(peak[0] == s2_.getPeak(i)[0]) {
					s2_.removePeak(i--);
				}
			}
		}
		
		ArrayList<double[]> experimentalPeaks = new ArrayList<double[]>();
		ArrayList<double[]> syntheticPeaks = new ArrayList<double[]>();
		
		int charge = s1_.getCharge();
		if(charge == 1) {
			charge = 2;
		}
		for(int c=1; c<charge; c++) {
			double[] peaks = p.getTheoreticalLadder(ProteomeConstants.Y_ION, c, false);
			for(int i=0; i<peaks.length; i++) {
				ArrayList<double[]> peakList = s1_.getSubPeaks(peaks[i]-0.02, peaks[i]+0.02);
				if(peakList.size() != 0) {
					double[] peak = {0,0};
					double closestMZ = 0;
					double intensity = 0;
					for(double[] thisPeak : peakList) {
						if(Math.abs(closestMZ - peaks[i]) > Math.abs(thisPeak[0] - peaks[i])) {
							closestMZ = peaks[i];
						}
						intensity += thisPeak[1];
					}
					peak[0] = closestMZ;
					peak[1] = intensity;
					experimentalPeaks.add(peak);
				}
				
				peakList = s2_.getSubPeaks(peaks[i]-0.02, peaks[i]+0.02);
				if(peakList.size() != 0) {
					double[] peak = {0,0};
					double closestMZ = 0;
					double intensity = 0;
					for(double[] thisPeak : peakList) {
						if(Math.abs(closestMZ - peaks[i]) > Math.abs(thisPeak[0] - peaks[i])) {
							closestMZ = peaks[i];
						}
						intensity += thisPeak[1];
					}
					peak[0] = closestMZ;
					peak[1] = intensity;
					syntheticPeaks.add(peak);
				}
			}
			
			peaks = p.getTheoreticalLadder(ProteomeConstants.B_ION, c, false);
			for(int i=0; i<peaks.length; i++) {
				ArrayList<double[]> peakList = s1_.getSubPeaks(peaks[i]-0.02, peaks[i]+0.02);
				if(peakList.size() != 0) {
					double[] peak = {0,0};
					double closestMZ = 0;
					double intensity = 0;
					for(double[] thisPeak : peakList) {
						if(Math.abs(closestMZ - peaks[i]) > Math.abs(thisPeak[0] - peaks[i])) {
							closestMZ = peaks[i];
						}
						intensity += thisPeak[1];
					}
					peak[0] = closestMZ;
					peak[1] = intensity;
					experimentalPeaks.add(peak);
				}
				
				peakList = s2_.getSubPeaks(peaks[i]-0.02, peaks[i]+0.02);
				if(peakList.size() != 0) {
					double[] peak = {0,0};
					double closestMZ = 0;
					double intensity = 0;
					for(double[] thisPeak : peakList) {
						if(Math.abs(closestMZ - peaks[i]) > Math.abs(thisPeak[0] - peaks[i])) {
							closestMZ = peaks[i];
						}
						intensity += thisPeak[1];
					}
					peak[0] = closestMZ;
					peak[1] = intensity;
					syntheticPeaks.add(peak);
				}
			}
		}
		
		// Remove duplicated peaks
		Hashtable<String, String> isDuplicatedPeak = new Hashtable<String, String>();
		for(int i=0; i<experimentalPeaks.size(); i++) {
			double[] peak = experimentalPeaks.get(i);
			if(isDuplicatedPeak.get(peak[1]+"") == null) {
				isDuplicatedPeak.put(peak[1]+"", "");
			} else {
				experimentalPeaks.remove(i--);
			}
 		}
		isDuplicatedPeak.clear();
		for(int i=0; i<syntheticPeaks.size(); i++) {
			double[] peak = syntheticPeaks.get(i);
			if(isDuplicatedPeak.get(peak[1]+"") == null) {
				isDuplicatedPeak.put(peak[1]+"", "");
			} else {
				syntheticPeaks.remove(i--);
			}
 		}
		isDuplicatedPeak.clear();
		
		if(experimentalPeaks.size() != 0 && syntheticPeaks.size() != 0) {
			Collections.sort(experimentalPeaks, new Comparator<double[]>() {
				@Override
				public int compare(double[] o1, double[] o2) {
					if(o1[0] < o2[0]) {
						return -1;
					}else if(o1[0] > o2[0]) {
						return 1;
					}
					return 0;
				}
			});
			Collections.sort(syntheticPeaks, new Comparator<double[]>() {
				@Override
				public int compare(double[] o1, double[] o2) {
					if(o1[0] < o2[0]) {
						return -1;
					}else if(o1[0] > o2[0]) {
						return 1;
					}
					return 0;
				}
			});
		} else {
			return 0;
		}
		
		double maxExpInt = 0;
		double maxThrInt = 0;
		for(int i=0; i<experimentalPeaks.size(); i++) {
//			maxExp = Math.max(experimentalPeaks.get(i)[1], maxExp);
			maxExp += experimentalPeaks.get(i)[1] * experimentalPeaks.get(i)[1];
			maxExpInt = Math.max(experimentalPeaks.get(i)[1], maxExpInt);
//			maxExp = Math.max(maxExp, s1_.getPeak(i)[1]);
//			maxExp += s1_.getPeak(i)[1];
		}
		for(int i=0; i<syntheticPeaks.size(); i++) {
//			maxThr = Math.max(syntheticPeaks.get(i)[1], maxThr);
			maxThr += syntheticPeaks.get(i)[1] * syntheticPeaks.get(i)[1];
			maxThrInt = Math.max(syntheticPeaks.get(i)[1], maxThrInt);
//			maxThr = Math.max(maxThr, s2_.getPeak(i)[1]);
//			maxThr += s2_.getPeak(i)[1];
		}

		for(int i=0; i<s1_.peaks.size(); i++) {
			if(s1_.peaks.get(i)[1] > maxExpInt) {
				s1_.peaks.remove(i--);
			}
		}
		
		for(int i=0; i<s2_.peaks.size(); i++) {
			if(s2_.peaks.get(i)[1] > maxThrInt) {
				s2_.peaks.remove(i--);
			}
		}
		
		//s1_.peaks = experimentalPeaks;
		//s2_.peaks = syntheticPeaks;
		
		int sizeOfExp = experimentalPeaks.size();
		int sizeOfThr = syntheticPeaks.size();
		double[][] expPeaks = new double[sizeOfExp][2];
		double[][] thrPeaks = new double[sizeOfThr][2];

		for(int i=0; i<sizeOfThr; i++) {
			thrPeaks[i][0] = syntheticPeaks.get(i)[0];
			thrPeaks[i][1] = syntheticPeaks.get(i)[1];
			thrPeaks[i][1] = thrPeaks[i][1] / Math.sqrt(maxThr);
		}
		
		for(int i=0; i<sizeOfExp; i++) {
			expPeaks[i][0] = experimentalPeaks.get(i)[0];
			expPeaks[i][1] = experimentalPeaks.get(i)[1];
			expPeaks[i][1] = expPeaks[i][1] / Math.sqrt(maxExp);
		}
		
		
		
		// normalization theoretical peaks.
		/*
		double norThr = 0;
		for(int i=0; i<sizeOfThr; i++) {
			thrPeaks[i][1] = thrPeaks[i][1];
			if(isSquareRootNorm) {
				thrPeaks[i][1] = Math.sqrt(thrPeaks[i][1]);
			}
			norThr += Math.pow(thrPeaks[i][1], 2);
		}
		norThr = Math.sqrt(norThr);
		for(int i=0; i<sizeOfThr; i++) thrPeaks[i][1] = thrPeaks[i][1] / norThr;
		
		// normalization experimental peaks.
		double norExp = 0;
		for(int i=0; i<sizeOfExp; i++) {
			expPeaks[i][1] = expPeaks[i][1];
			if(isSquareRootNorm) {
				expPeaks[i][1] = Math.sqrt(expPeaks[i][1]);
			}
			norExp += Math.pow(expPeaks[i][1],2);
		}
		norExp = Math.sqrt(norExp);
		for(int i=0; i<sizeOfExp; i++) expPeaks[i][1] = expPeaks[i][1] / norExp;
		*/
		// 
		double innerproduct = 0;
		// inner product
		int thrIndex = 0;
		int expIndex = 0;
		boolean doMatchFurther = true;
		while(doMatchFurther) {
			double delta = thrPeaks[thrIndex][0] - expPeaks[expIndex][0];
			if(Math.abs(delta) <= tolerance) {
				innerproduct += thrPeaks[thrIndex][1] * expPeaks[expIndex][1];
				expIndex++;
				thrIndex++;
			}
			else if(delta > 0) expIndex++;
			else if(delta < 0) thrIndex++;
			else if(delta == 0) {
				expIndex++;
				thrIndex++;
			}
			
			if(thrIndex >= thrPeaks.length || expIndex >= expPeaks.length) doMatchFurther = false;
		}
		
		// calculate experimental Norm.
		DecimalFormat decimalFormat = new DecimalFormat("#.#####");
		innerproduct = Double.parseDouble(decimalFormat.format(innerproduct));
		double sca = 1 - 2 * (Math.acos(innerproduct) / Math.PI);
		
		return sca;
	}
	
	public static double spectralCorrelation (ArrayList<double[]> experimentalPeaks, ArrayList<double[]> theoreticalPeaks, double tolerance) {
		int sizeOfExp = experimentalPeaks.size();
		int sizeOfThr = theoreticalPeaks.size();
		double[][] expPeaks = new double[sizeOfExp][2];
		double[][] thrPeaks = new double[sizeOfThr][2];
		
		// normalization theoretical peaks.
		double norThr = 0;
		for(int i=0; i<sizeOfThr; i++) {
			thrPeaks[i][0] = theoreticalPeaks.get(i)[0];
			thrPeaks[i][1] = Math.sqrt(theoreticalPeaks.get(i)[1]);
			norThr += theoreticalPeaks.get(i)[1];
		}
		norThr = Math.sqrt(norThr);
		for(int i=0; i<sizeOfThr; i++) thrPeaks[i][1] = thrPeaks[i][1] / norThr;
		
		// normalization experimental peaks.
		double norExp = 0;
		for(int i=0; i<sizeOfExp; i++) {
			expPeaks[i][0] = experimentalPeaks.get(i)[0];
			expPeaks[i][1] = Math.sqrt(experimentalPeaks.get(i)[1]);
			norExp += experimentalPeaks.get(i)[1];
		}
		norExp = Math.sqrt(norExp);
		for(int i=0; i<sizeOfExp; i++) expPeaks[i][1] = expPeaks[i][1] / norExp;
		
		// 
		double innerproduct = 0;
		// inner product
		int thrIndex = 0;
		int expIndex = 0;
		boolean doMatchFurther = true;
		while(doMatchFurther) {
			double delta = thrPeaks[thrIndex][0] - expPeaks[expIndex][0];
			if(Math.abs(delta) <= tolerance) {
				innerproduct += thrPeaks[thrIndex][1] * expPeaks[expIndex][1];
				expIndex++;
				thrIndex++;
			}
			else if(delta > 0) expIndex++;
			else if(delta < 0) thrIndex++;
			else if(delta == 0) {
				expIndex++;
				thrIndex++;
			}
			
			if(thrIndex >= thrPeaks.length || expIndex >= expPeaks.length) doMatchFurther = false;
		}
		
		// calculate experimental Norm.
		DecimalFormat decimalFormat = new DecimalFormat("#.####");
		innerproduct = Double.parseDouble(decimalFormat.format(innerproduct));
		double sca = 1 - 2 * (Math.acos(innerproduct) / Math.PI);
		
		return sca;
	}
}