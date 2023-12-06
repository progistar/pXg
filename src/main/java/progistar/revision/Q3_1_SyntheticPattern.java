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
import progistar.thirdparty.sSim.Spectrum;

class Target_ {
	String title;
	int charge;
	String fullRecord;
	Spectrum spectrum;
}

public class Q3_1_SyntheticPattern {

	public static void main(String[] args) throws IOException {
		File file = new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/5.Synthetic/synthetic_list.txt");
		File[] files = new File("/Users/gistar/projects/pXg/MSMS").listFiles();
		Hashtable<String, ArrayList<Target_>> targets = getSyntheticList(file, files);	
		
		Iterator<String> peptides = (Iterator<String>) targets.keys();
		
		while(peptides.hasNext()) {
			String peptide = peptides.next();
			BufferedWriter BW = new BufferedWriter(new FileWriter("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/"+peptide+"_exp.mgf"));
			ArrayList<Target_> list = targets.get(peptide);
			
			for(Target_ t : list) {
				BW.append("BEGIN IONS");
				BW.newLine();
				
				BW.append("TITLE=").append(t.spectrum.getTitle());
				BW.newLine();
				BW.append("RTINSECONDS=").append(t.spectrum.getRT()+"");
				BW.newLine();
				BW.append("PEPMASS=").append(t.spectrum.getPrecursorMz()+" "+t.spectrum.getPrecursorInt());
				BW.newLine();
				BW.append("CHARGE=").append(t.spectrum.getCharge()+"");
				BW.newLine();
				
				for(double[] peak : t.spectrum.peaks) {
					BW.append(peak[0]+" "+peak[1]);
					BW.newLine();
				}
				BW.append("END IONS");
				BW.newLine();
			}
			
			BW.close();
		}
		
		
	}
	
	public static Hashtable<String, ArrayList<Target_>> getSyntheticList (File file, File[] mgfs) throws IOException {
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
		
		Hashtable<String, ArrayList<Target_>> targets = new Hashtable<String, ArrayList<Target_>>();	
		
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
			
			ArrayList<Target_> targets_ = targets.get(peptide);
			Target_ target = new Target_();
			target.charge = Integer.parseInt(charge);
			target.fullRecord = line;
			String title = specId.split("\\.")[0]+"."+scanId+"."+scanId+"."+charge;
			target.title = title;
			if(targets_ == null) {
				targets_ = new ArrayList<Target_>();
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
