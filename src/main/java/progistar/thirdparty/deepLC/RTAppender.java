package progistar.thirdparty.deepLC;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class RTAppender {

	public static Pattern MODREG = Pattern.compile("([+0-9.]+)");
	public static String MET = "Oxidation";
	public static String DEAMI = "Deamidated";
	public static String CARBAM = "Carbamidomethyl";
	public static String PHOSPHO = "Phospho";
	public static String CYSTEINLY = "Cysteinyl";
	
	public static void main(String[] args) throws IOException {
		String pXgFileName = "C:\\Users\\progi\\Desktop\\Projects\\pXg\\Qi_MCP2021\\Results\\h1975wc.pxg.netMHCpan";
		String deeplcFolder = "C:\\Users\\progi\\Desktop\\Projects\\pXg\\Qi_MCP2021\\RT";
		
		File[] files = new File(deeplcFolder).listFiles();
		// key
		// fileName_InferredPeptide_PTM_RT
		// value
		// predicted RT
		Hashtable<String, String> RTMapper = new Hashtable<String, String>();
		
		for(File file : files) {
			if(file.getName().contains("deeplc_predictions")) {
				BufferedReader BR = new BufferedReader(new FileReader(file));
				String line = null;
				
				String fileName = file.getName().split("\\.")[0]+".raw";
				BR.readLine(); // skip header
				while((line = BR.readLine()) != null) {
					String[] fields = line.split("\\,");
					String key = fileName+"_"+fields[1]+"_"+fields[2]+"_"+Double.parseDouble(fields[3]);
					
					RTMapper.put(key, fields[4]);
				}
				BR.close();
			}
		}
		
		// read PXG	
		Hashtable<String, String> massToModName = new Hashtable<String, String>();
		massToModName.put("+15.99", MET);
		massToModName.put("+.98", DEAMI);
		massToModName.put("+119.00", CYSTEINLY);
		massToModName.put("+79.97", PHOSPHO);
		massToModName.put("+58.01", CARBAM);
		
		int peaksPeptideIndex = 3;
		int rtIndex = 11;
		int fractionIndex = 1;
		int peptideIndex = 19;
		
		BufferedReader BR = new BufferedReader(new FileReader(pXgFileName));
		
		String line = null;
		
		// header
		System.out.println(BR.readLine()+"\tdeeplcRT");
		
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String modifications = "";
			if(fields[peaksPeptideIndex].contains("+") || fields[peaksPeptideIndex].contains("-")) {
				Matcher matcher = MODREG.matcher(fields[peaksPeptideIndex]);
				while(matcher.find()) {
					String group = matcher.group();
					int pos = fields[peaksPeptideIndex].indexOf("("+group+")");
					fields[peaksPeptideIndex] = fields[peaksPeptideIndex].replaceFirst("\\([+0-9.]+\\)", "");
					modifications += pos+"|"+massToModName.get(group)+"|";
				}
				modifications = modifications.substring(0, modifications.length()-1);
			}
			//
			
			String key = fields[fractionIndex]+"_"+fields[peptideIndex]+"_"+modifications+"_"+Double.parseDouble(fields[rtIndex]);
			String pRT = RTMapper.get(key);
			
			if(pRT == null) {
				// multiple genomic loci with ambiguous I/L
				System.out.println(line+"\t-");
			} else {
				System.out.println(line+"\t"+pRT);
			}
		}
		
		
		BR.close();
	}
}
