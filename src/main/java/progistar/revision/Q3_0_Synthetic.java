package progistar.revision;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;

import progistar.thirdparty.sSim.Spectra;
import progistar.thirdparty.sSim.SpectralScores;
import progistar.thirdparty.sSim.Spectrum;

class Target {
	String title;
	int charge;
	String fullRecord;
	Spectrum spectrum;
}

public class Q3_0_Synthetic {

	public static void main(String[] args) throws IOException {
		File file = new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/5.Synthetic/synthetic_list.txt");
		File[] files = new File("/Users/gistar/projects/pXg/MSMS").listFiles();
		File[] syntheticFiles = new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/5.Synthetic/Synthetic_MGF").listFiles();
		Hashtable<String, ArrayList<Target>> targets = getSyntheticList(file, files);	
		Hashtable<String, Spectra> syntheticSpectra = new Hashtable<String, Spectra>();
		
		for(File syntheticFile : syntheticFiles) {
			String peptide = syntheticFile.getName().split("\\_")[0];
			Spectra spectra = new Spectra(syntheticFile.getAbsolutePath(), Spectra.FILE_TYPE_MGF);
			
			syntheticSpectra.put(peptide, spectra);
		}
		
		Iterator<String> peptides = (Iterator<String>) targets.keys();
		BufferedWriter BW = new BufferedWriter(new FileWriter("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/5.Synthetic/CompareResults.tsv"));
		while(peptides.hasNext()) {
			String peptide = peptides.next();
			ArrayList<Target> targets_ = targets.get(peptide);
			Spectra spectra = syntheticSpectra.get(peptide);
			
			double sca = 0;
			double fullSCA = 0;
			double pcc = 0;
			Spectrum bestExp = null;
			Spectrum bestSyn = null;
			for(Target target : targets_) {
				int size = spectra.sizeOfSpectra();
				for(int i=0; i<size;i ++) {
					Spectrum spectrum = spectra.getSpectrumByIndex(i);
					spectrum.setPeptide(peptide);
					target.spectrum.setPeptide(peptide);
					
					if(spectrum.getCharge() == 0) continue;
					if(spectrum.getCharge() != target.spectrum.getCharge()) continue;
					
					
					double thisPCC = SpectralScores.spectralCorrelation(target.spectrum.deepCopy(), spectrum.deepCopy(), 0.02);
					double thisSCA = SpectralScores.spectralContrastAngle(target.spectrum.deepCopy(), spectrum.deepCopy(), 0.02, false);
					double thisFullSCA = SpectralScores.spectralContrastAngleWithFull(target.spectrum.deepCopy(), spectrum.deepCopy(), 0.02, false);
					if((sca + fullSCA + pcc) < (thisPCC + thisSCA + thisFullSCA)) {
						pcc = thisPCC;
						sca = thisSCA;
						fullSCA = thisFullSCA;
						
						bestExp = target.spectrum;
						bestSyn = spectrum;
					}
				}
			}
			
			Spectrum[] two = {bestExp, bestSyn};
			BW.append("PEPTIDE: "+peptide);
			BW.newLine();
			BW.append("PCC: "+pcc+" SA: "+sca+" FullSA: "+fullSCA);
			BW.newLine();
			for(Spectrum spec : two ) {
				BW.append(spec.getTitle());
				BW.newLine();
				
				for(double[] peak : spec.peaks) {
					BW.append(peak[0]+" "+peak[1]);
					BW.newLine();
				}
			}
			System.out.println(peptide+":"+pcc);
			System.out.println(peptide+":"+sca);
			System.out.println(peptide+":full "+fullSCA);
		}
		BW.close();
		
		
		
	}
	
	public static Hashtable<String, ArrayList<Target>> getSyntheticList (File file, File[] mgfs) throws IOException {
		ArrayList<Spectra> specSets = new ArrayList<Spectra>();
		for(int i = 0; i < mgfs.length; i++) {
			if(mgfs[i].getName().startsWith("THP1_2") || 
					mgfs[i].getName().startsWith("THP1_3") ||
					mgfs[i].getName().startsWith("B_LCL1") ||
					mgfs[i].getName().startsWith("B_LCL2") ||
					mgfs[i].getName().startsWith("B_LCL3") ||
					mgfs[i].getName().startsWith("B_LCL4"))
			System.out.println(" read file: "+mgfs[i].getAbsolutePath());
			specSets.add(new Spectra(mgfs[i].getAbsolutePath(), Spectra.FILE_TYPE_MGF));
		}
		
		Hashtable<String, ArrayList<Target>> targets = new Hashtable<String, ArrayList<Target>>();	
		
		// read experimental peptides
		int specIdIdx = 1;
		int scanIdx = 6;
		int chargeIdx = 12;
		int peptideIdx = 23;
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = BR.readLine(); // skip header
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[peptideIdx];
			String specId = fields[specIdIdx];
			String charge = fields[chargeIdx];
			String scanId = fields[scanIdx].split("\\:")[1];
			
			ArrayList<Target> targets_ = targets.get(peptide);
			Target target = new Target();
			target.charge = Integer.parseInt(charge);
			target.fullRecord = line;
			String title = specId.split("\\.")[0]+"."+scanId+"."+scanId+"."+charge;
			target.title = title;
			if(targets_ == null) {
				targets_ = new ArrayList<Target>();
				targets.put(peptide, targets_);
			}
			targets_.add(target);
			Spectrum spectrum = null;
			for(Spectra spectra : specSets) {
				spectrum = spectra.getSpectrumByScanNum(title);
				if(spectrum != null) {
					break;
				}
			}
			if(spectrum == null) {
				System.out.println("We cannot find: "+title);
			} else {
				target.spectrum = spectrum;
			}
 		}
		
		BR.close();
		
		System.out.println(targets.size());
		
		return targets;
		
	}
}
