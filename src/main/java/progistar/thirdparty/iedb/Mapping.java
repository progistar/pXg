package progistar.thirdparty.iedb;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

public class Mapping {

	public static void main(String[] args) throws IOException {
		File previousReportFile = new File("/Users/gistar/projects/pXg/PreviousStudyResource/IEDB_IEAtlas_HLAligand_Cuevas_Laumont_Scull.txt");
		File pXgReportFile = new File("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/canonical_fdr5.tsv");
		
		Hashtable<String, String> mapper = loadPreviousMAPList(previousReportFile);
		String[] sources = {
				"IEDB",
				"IEDB T Cell+",
				"IEDB T Cell-",
				"IEAtlas cancer",
				"IEAtlas normal",
				"HLA-ligand",
				"Canonical (Scull et al.)",
				"Canonical (Cuevas et al.)",
				"Canonical (Laumont et al.)",
				"Noncanonical (Scull et al.)",
				"Noncanonical (Cuevas et al.)",
				"Noncanonical (Laumont et al.)"
		};
		
		BufferedReader BR = new BufferedReader(new FileReader(pXgReportFile));
		BufferedWriter BW = new BufferedWriter(new FileWriter(pXgReportFile.getAbsolutePath().replace(".tsv", ".reported.tmp.tsv")));
		String line = null;
		
		String header = BR.readLine();
		BW.append(header);
		for(int i=0; i<sources.length; i++) {
			BW.append("\t"+sources[i]);
		}
		BW.append("\tNumberOfMappedSources");
		BW.newLine();
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[24];
			
			String record = line;
			int mappingCount= 0;
			for(int i=0; i<sources.length; i++) {
				String key = peptide+"_"+sources[i];
				if(mapper.get(key) == null) {
					record +="\tNo";
				} else {
					record +="\tYes";
					mappingCount++;
				}
			}
			record+="\t"+mappingCount;
			BW.append(record);
			BW.newLine();
		}
		
		BW.close();
		BR.close();
		
	}
	
	public static Hashtable<String, String> loadPreviousMAPList (File file) throws IOException {
		Hashtable<String, String> mapper = new Hashtable<String, String>();
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		Hashtable<String, String> IEDB = new Hashtable<String, String>();
		Hashtable<String, String> HLAligand = new Hashtable<String, String>();
		Hashtable<String, String> Scull = new Hashtable<String, String>();
		Hashtable<String, String> Cuevas = new Hashtable<String, String>();
		Hashtable<String, String> Laumont = new Hashtable<String, String>();
		Hashtable<String, String> IEAtlas = new Hashtable<String, String>();
		Hashtable<String, String> BigDBs = new Hashtable<String, String>();
		Hashtable<String, String> ThreeDBs = new Hashtable<String, String>();
		Hashtable<String, String> AllDBs = new Hashtable<String, String>();
		
		BR.readLine(); // skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[0];
			String source = fields[2];
			
			mapper.put(peptide+"_"+source, "");
			
			if(source.contains("IEDB")) {
				IEDB.put(peptide, "");
			}else if(source.contains("IEAtlas")) {
				IEAtlas.put(peptide, "");
			}else if(source.contains("Scull")) {
				Scull.put(peptide, "");
			}else if(source.contains("Laumont")) {
				Laumont.put(peptide, "");
			}else if(source.contains("Cuevas")) {
				Cuevas.put(peptide, "");
			}else if(source.contains("HLA-ligand")) {
				HLAligand.put(peptide, "");
			}
		}
		
		BR.close();
		BigDBs.putAll(IEAtlas);
		BigDBs.putAll(IEDB);
		BigDBs.putAll(HLAligand);
		System.out.println("IEDB: "+IEDB.size());
		System.out.println("HLA-ligand: "+HLAligand.size());
		System.out.println("IEAtlas: "+IEAtlas.size());
		System.out.println("Union above: "+BigDBs.size());
		
		
		System.out.println("Scull: "+Scull.size());
		System.out.println("Cuevas: "+Cuevas.size());
		System.out.println("Laumont: "+Laumont.size());
		ThreeDBs.putAll(Scull);
		ThreeDBs.putAll(Cuevas);
		ThreeDBs.putAll(Laumont);
		AllDBs.putAll(BigDBs);
		AllDBs.putAll(ThreeDBs);
		System.out.println("Union publications: "+ThreeDBs.size());
		System.out.println("Union: "+AllDBs.size());
		
		return mapper;
	}
}
