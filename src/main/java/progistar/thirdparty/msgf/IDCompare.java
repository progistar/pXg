package progistar.thirdparty.msgf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

public class IDCompare {

	public static final int MSGF_SPEC_FILE_IDX = 0;
	public static final int MSGF_SCAN_IDX = 2;
	public static final int MSGF_PEPTIDE_IDX = 9;
	
	public static final int PXG_SPEC_FILE_IDX = 1;
	public static final int PXG_SCAN_IDX = 4;
	public static final int PXG_PEPTIDE_IDX = 20;
	public static final int PXG_IS_CANONICAL = 36;
	
	public static Hashtable<String, String> readPXG (File file) throws IOException {
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		BR.readLine(); // skip header
		
		Hashtable<String, String> pxgPSMs = new Hashtable<String, String>();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String key = fields[PXG_SPEC_FILE_IDX] +"_"+fields[PXG_SCAN_IDX].split("\\:")[1];
			String peptide = fields[PXG_PEPTIDE_IDX];
			peptide = peptide.replaceAll("[+-.\\(\\)0123456789*]", "");
			
			pxgPSMs.put(key, fields[3]+"_"+peptide);
		}
		
		BR.close();
		
		return pxgPSMs;
	}
	
	public static Hashtable<String, String> readMSGF (File file) throws IOException {
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		Hashtable<String, String> msgfPSMs = new Hashtable<String, String>();
		
		BR.readLine(); // skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String key = fields[MSGF_SPEC_FILE_IDX] +"_"+fields[MSGF_SCAN_IDX];
			String peptide = fields[MSGF_PEPTIDE_IDX];
			peptide = peptide.replaceAll("[+-.0123456789*]", "");
			msgfPSMs.put(key, fields[MSGF_PEPTIDE_IDX]+"_"+peptide);
		}
		
		BR.close();
		
		return msgfPSMs;
	}
	
	public static void main(String[] args) throws IOException {
		File msgfRes = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/MSGF/fdr_at5/S3.fdr");
		File pxgRes = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/3.withCalibrationAddScanNum/pXg/PeptideAnnotationS3_5ppm_002_recal.scanNum.rep1.rank10.pXg.BA.fdr");
		File pxgMetOnlyRes = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/5.withCalibrationAddScanNumWithoutDeami/pXg/PeptideAnnotationS3_5ppm_002_recal.scanNum.rep1.rank10.pXg.fdr");
		
		Hashtable<String, String> msgfPSMs = readMSGF(msgfRes);
		Hashtable<String, String> pxgPSMs = readPXG(pxgRes);
		Hashtable<String, String> pxgMetOnlyPSMs = readPXG(pxgMetOnlyRes);
		
		// 0: all same
		// 1: pxg same
		// 2: pxgMet same
		// 3: not same but pxg and pxgMet are same
		// 4: totally different among them
		System.out.println("Scan,MSGFPept,PEAKSwDeami,PEAKSwoDeami,ClassCode");
		
		msgfPSMs.forEach((scan, pept)->{
			String fullPept = pept.split("_")[0];
			String stripPept = pept.split("_")[1];
			
			String pxgPept = pxgPSMs.get(scan);
			String fullPXGPept = "NA";
			String stripPXGPept = "NA";
			
			String pxgMetPept = pxgMetOnlyPSMs.get(scan);
			String fullPXGMetPept = "NA";
			String stripPXGMetPept = "NA";
			
			boolean isPXGSame = false;
			boolean isPXGMetSame = false;
			boolean isPXGMetPXGSame = false;
			
			if(pxgPept != null) {
				fullPXGPept = pxgPept.split("_")[0];
				stripPXGPept = pxgPept.split("_")[1];
				
				if(stripPXGPept.equalsIgnoreCase(stripPept)) {
					isPXGSame = true;
				}
			}
			
			if(pxgMetPept != null) {
				fullPXGMetPept = pxgMetPept.split("_")[0];
				stripPXGMetPept = pxgMetPept.split("_")[1];
				
				if(stripPXGMetPept.equalsIgnoreCase(stripPept)) {
					isPXGMetSame = true;
				}
			}
			
			if(pxgPept != null && pxgMetPept != null) {
				if(stripPXGPept.equalsIgnoreCase(stripPXGMetPept)) {
					isPXGMetPXGSame = true;
				}
			}
			
			int idx = 0;
			// all same
			if(isPXGSame && isPXGMetSame) {
				idx = 0;
			} else if(isPXGSame && !isPXGMetSame) {
				idx = 1;
			} else if(!isPXGSame && isPXGMetSame) {
				idx = 2;
			} else if(!isPXGSame && !isPXGMetSame && isPXGMetPXGSame) {
				idx = 3;
			} else {
				idx = 4;
			}
			
			System.out.println(scan+","+fullPept+","+fullPXGPept+","+fullPXGMetPept+","+idx);
			
		});
	}
}
