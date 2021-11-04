package progistar.thirdparty.deepLC;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class pXg2DeepLCInput {

	public static Pattern MODREG = Pattern.compile("([+0-9.]+)");
	public static String MET = "Oxidation";
	public static String DEAMI = "Deamidated";
	public static String CARBAM = "Carbamidomethyl";
	public static String PHOSPHO = "Phospho";
	public static String CYSTEINLY = "Cysteinyl";
	
	public static void main(String[] args) throws IOException {
		
		Hashtable<String, String> massToModName = new Hashtable<String, String>();
		massToModName.put("+15.99", MET);
		massToModName.put("+.98", DEAMI);
		massToModName.put("+119.00", CYSTEINLY);
		massToModName.put("+79.97", PHOSPHO);
		massToModName.put("+58.01", CARBAM);
		
		int peaksPeptideIndex = 3;
		int rtIndex = 11;
		int fractionIndex = 1;
		int scanIndex = 4;
		int peptideIndex = 19;
		
		String fileName = "C:\\Users\\progi\\Desktop\\Projects\\pXg\\Laumont_NatCommun2016\\Results\\subjectM.5ppm.002.rep1.pXg";
		
		BufferedReader BR = new BufferedReader(new FileReader(fileName));
		
		String line = null;
		
		
		BR.readLine(); // skip header
		Hashtable<String, String> rmDuplications = new Hashtable<String, String>();
		Hashtable<String, ArrayList<String>> fractions = new Hashtable<String, ArrayList<String>>();
		
		
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
			
			ArrayList<String> fraction = fractions.get(fields[fractionIndex]);
			if(fraction == null) {
				fraction = new ArrayList<String>();
				fractions.put(fields[fractionIndex], fraction);
			}
			
			if(rmDuplications.get(fields[fractionIndex]+"_"+fields[scanIndex]+"_"+fields[peptideIndex]) == null) {
				
				String key = fields[peptideIndex]+","+modifications+","+fields[rtIndex];
				fraction.add(key);
				rmDuplications.put(fields[fractionIndex]+"_"+fields[scanIndex]+"_"+fields[peptideIndex], "");
			}
		}
		
		fractions.forEach( (fName, fraction) -> {
			try {
				BufferedWriter BW = new BufferedWriter(new FileWriter(fName+".csv"));
				BW.append("seq,modifications,tr");
				BW.newLine();
				
				for(String id : fraction) {
					BW.append(id);
					BW.newLine();
				}
				
				BW.close();
			}catch(IOException ioe) {
				
			}
		});
		
		
		BR.close();
	}
}
