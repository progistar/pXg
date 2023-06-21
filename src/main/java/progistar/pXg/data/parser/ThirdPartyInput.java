package progistar.pXg.data.parser;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;

public class ThirdPartyInput {
	
	public static String CE = null;
	
	public static String pxgResultPath = null;
	public static int rtIdx = -1;
	public static int infPeptIdx = -1;
	public static int isCanonicalIdx = -1;
	public static int labelIdx = -1;
	public static int chargeIdx = -1;
	public static int scoreIdx = -1;

	public static void printUsage() {
		System.out.println("Usage");
		System.out.println();
		System.out.println("Mandatory Fields");
		System.out.println("  -pxg                 : pXg result file path. We recommand to use the same gtf corresponding to alignment.");
		System.out.println("  -rt_col              : Retention time column index in the pXg result. One-based!");
		System.out.println("  -ipept_col           : Inffererd peptide column index in the pXg result. One-based!");
		System.out.println("  -ic_col              : IsCanonical column index in the pXg result. One-based!");
		System.out.println("  -label_col           : Label column index in the pXg result. One-based!");
		System.out.println("  -charge_col          : Charge column index in the pXg result. One-based!");
		System.out.println("  -score_col           : Score column index in the pXg result. One-based!");
		System.out.println("  -ce                  : SAM file path. The sam file must be sorted by coordinate.");
		
		System.out.println("Example1");
		System.out.println("java -Xmx30G -c pXg.jar progistar.pXg.data.parser.ThridPartyInput -pxg result.pxg -rt_col 13 -ipept_col 23 -ic_col 39 -label_col 1 -charge_col 12 -score_col 9 -ce 35");
	}
	
	public static boolean checkInput(String[] args) {
		// input check
		System.out.println(Constants.VERSION+" "+Constants.RELEASE);
		System.out.println(Constants.INTRODUCE);
		System.out.println();
		

		
		// print parameter description
		if(args.length == 0) {
			printUsage();
			return false;
		}
		
		for(int i=0; i<args.length; i+=2) {
			String option = args[i].toLowerCase();
			
			// -gtf (mandatory)
			if(option.equalsIgnoreCase("-pxg")) {
				pxgResultPath = args[i+1];
				
			} else if(option.equalsIgnoreCase("-rt_col")) {
				rtIdx = Integer.parseInt(args[i+1]);
				
			} else if(option.equalsIgnoreCase("-ipept_col")) {
				infPeptIdx = Integer.parseInt(args[i+1]);
				
			} else if(option.equalsIgnoreCase("-ic_col")) {
				isCanonicalIdx = Integer.parseInt(args[i+1]);
				
			} else if(option.equalsIgnoreCase("-label_col")) {
				labelIdx = Integer.parseInt(args[i+1]);
				
			} else if(option.equalsIgnoreCase("-charge_col")) {
				chargeIdx = Integer.parseInt(args[i+1]);
				
			} else if(option.equalsIgnoreCase("-score_col")) {
				scoreIdx = Integer.parseInt(args[i+1]);
				
			} else if(option.equalsIgnoreCase("-ce")) {
				CE = args[i+1];
			}
		}
		
		if(pxgResultPath != null &&
				rtIdx != -1 &&
				infPeptIdx != -1 &&
				isCanonicalIdx != -1 &&
				labelIdx != -1 &&
				chargeIdx != -1 &&
				scoreIdx != -1 &&
				CE != null) {
			return true;
		}
		printUsage();
		return false;
	}
	
	public static void main(String[] args) throws IOException {
		/*
		if(!checkInput(args)) {
			return;
		}
		*/
		
		boolean[] selectedSet = {
				false, // B-LCL1
				false, // B-LCL2
				false, // B-LCL3
				false, // B-LCL4
				false, // DOHH2
				false, // HBL1
				false, // SUDHL4
				false, // THP1-1
				false, // THP1-2
				false, // THP1-3
				
				true, // DI2T
				true, // DI5T
				true, // IN19T
				true, // IN26T
				true, // IN81T
				true, // R1_IN403T
				true, // R2_IN403T
				true, // IN407T
				true, // IN506T
				true, // IN524T
				true, // IN525T
				true, // IN526T
				true, // IN529T
				true, // R1_IN532T
				true, // R2_IN532T
				true, // M004T
				true // M009T
		};
		
		String[] pxgList = {
				"B-LCL1.pxg",
				"B-LCL2.pxg",
				"B-LCL3.pxg",
				"B-LCL4.pxg",
				"DOHH2.pxg",
				"HBL1.pxg",
				"SUDHL4.pxg",
				"THP1-1.pxg",
				"THP1-2.pxg",
				"THP1-3.pxg",
				
				"DI2T.pxg",
				"DI5T.pxg",
				"IN19T.pxg",
				"IN26T.pxg",
				"IN81T.pxg",
				"R1_IN403T.pxg",
				"R2_IN403T.pxg",
				"IN407T.pxg",
				"IN506T.pxg",
				"IN524T.pxg",
				"IN525T.pxg",
				"IN526T.pxg",
				"IN529T.pxg",
				"R1_IN532T.pxg",
				"R2_IN532T.pxg",
				"M004T.pxg",
				"M009T.pxg"
		};
		
		String absPath = "/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut";
		absPath = "/Users/gistar/projects/GastricCancer_NCC/pXg";
		// CE = 35 for B-LCL
		// CE = 25 for SHUDL4, DOHH2, HBL1
		// CE = 32 for THP
		
		String[] CESet = {
			"35",
			"35",
			"35",
			"35",
			"25",
			"25",
			"25",
			"32",
			"32",
			"32",
			
			"28", // DI2T
			"28", // DI5T
			"28", // IN19T
			"28", // IN26T
			"28", // IN81T
			"28", // R1_IN403T
			"28", // R2_IN403T
			"28", // IN407T
			"28", // IN506T
			"28", // R1_IN523T
			"28", // R2_IN523T
			"28", // IN524T
			"28", // IN525T
			"28", // IN526T
			"28", // IN529T
			"28", // M004T
			"28" // M009T
		};
		
		int rtIdx = 14;
		int infPeptIdx = 24;
		int isCanonicalIdx = 41;
		int labelIdx = 2;
		int chargeIdx = 13;
		int alcIdx = 10;
		
		String prositHeader = "modified_sequence,collision_energy,precursor_charge";
		String deepLCHeader = "seq,modifications,tr";
		for(int idx = 0; idx < pxgList.length; idx++) {
			if(!selectedSet[idx]) {
				continue;
			}
			
			String CE = CESet[idx];
			String pxg = pxgList[idx];
			File file = new File(absPath+"/"+pxg);
			BufferedReader BR = new BufferedReader(new FileReader(file));
			BufferedWriter BWProsit = new BufferedWriter(new FileWriter(file.getAbsolutePath().replace(".pxg", ".input.prosit")));
			BufferedWriter BWDeepLC = new BufferedWriter(new FileWriter(file.getAbsolutePath().replace(".pxg", ".input.deeplc")));
			BufferedWriter BWDeepLCCal = new BufferedWriter(new FileWriter(file.getAbsolutePath().replace(".pxg", ".input.cal.deeplc")));
			BufferedWriter BWNetMHCpan = new BufferedWriter(new FileWriter(file.getAbsolutePath().replace(".pxg", ".input.netmhcpan")));
			String line = null;
			
			Hashtable<String, String> prositInput = new Hashtable<String, String>();
			Hashtable<String, String> deepLCInput = new Hashtable<String, String>();
			Hashtable<String, String> deepLCInputCal = new Hashtable<String, String>();
			Hashtable<String, Double> isBestScore = new Hashtable<String, Double>();
			Hashtable<String, String> netmhcpanInput = new Hashtable<String, String>();
			BR.readLine();// skip header
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String prositList = fields[infPeptIdx]+","+CE+","+fields[chargeIdx];
				prositInput.put(prositList, "");
				String deepLCList = fields[infPeptIdx]+",,"+fields[rtIdx];
				deepLCInput.put(deepLCList,"");
				// canonical target sequence is input for the calibration
				// list only the best score
				if(fields[labelIdx].equalsIgnoreCase("1") && fields[isCanonicalIdx].equalsIgnoreCase("TRUE")) {
					
					Double alcScore = Double.parseDouble(fields[alcIdx]);
					if(isBestScore.get(fields[infPeptIdx]) == null) {
						deepLCInputCal.put(deepLCList,"");
						isBestScore.put(fields[infPeptIdx], alcScore);
					} else if(isBestScore.get(fields[infPeptIdx]) < alcScore) {
						deepLCInputCal.put(deepLCList,"");
						isBestScore.put(fields[infPeptIdx], alcScore);
					}
					
				}
				
				if(netmhcpanInput.get(fields[infPeptIdx]) == null) {
					BWNetMHCpan.append(fields[infPeptIdx]);
					BWNetMHCpan.newLine();
					
					netmhcpanInput.put(fields[infPeptIdx], "");
				}
			}
			
			BWProsit.append(prositHeader);
			BWProsit.newLine();
			prositInput.forEach((record, nil)->{
				try {
					BWProsit.append(record);
					BWProsit.newLine();
				}catch(IOException ioe) {
					
				}
			});
			
			BWDeepLC.append(deepLCHeader);
			BWDeepLC.newLine();
			deepLCInput.forEach((record, nil)->{
				try {
					BWDeepLC.append(record);
					BWDeepLC.newLine();
				}catch(IOException ioe) {
					
				}
			});
			
			BWDeepLCCal.append(deepLCHeader);
			BWDeepLCCal.newLine();
			deepLCInputCal.forEach((record, nil)->{
				try {
					BWDeepLCCal.append(record);
					BWDeepLCCal.newLine();
				}catch(IOException ioe) {
					
				}
			});
			
			
			BWProsit.close();
			BWDeepLC.close();
			BWDeepLCCal.close();
			BWNetMHCpan.close();
			BR.close();
		}
		
	}
}
