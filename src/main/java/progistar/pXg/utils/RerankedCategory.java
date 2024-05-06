package progistar.pXg.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

class PXGResult {
	public String record;
	public String peptide;
	public String topPeptide;
	public String isCanonical;
	public int rank;

	@Override
	public String toString () {
		int idx = getCategory();
		String category = "";
		if(idx == 0) {
			category = "top-ranked";
		} else if(idx == 1) {
			category = "scrambled";
		}  else if(idx == 2) {
			category = "deamidated";
		} else if(idx == 3) {
			category = "others";
		}

		return rank+"\t"+category+"\t"+isCanonical;
	}

	public int getCategory () {
		int category = 0; // same

		if(rank == 1) {
			return category; // same
		}

		int[] thisAAs = new int[26];
		int[] topAAs = new int[26];

		String thisStrip = peptide.replaceAll("[+-.0123456789\\(\\)*]*", "");
		String thisTopStrip = topPeptide.replaceAll("[+-.0123456789\\(\\)*]*", "");

		if(thisStrip.length() != thisTopStrip.length()) {
			return 3; // others
		}

		for(int i=0; i<thisStrip.length(); i++) {
			thisAAs[thisStrip.charAt(i)-'A']++;
		}

		for(int i=0; i<thisTopStrip.length(); i++) {
			topAAs[thisTopStrip.charAt(i)-'A']++;
		}

		boolean isScrambled = true;
		for(int i=0; i<thisStrip.length(); i++) {
			if(thisAAs[i] != topAAs[i]) {
				isScrambled = false;
			}
		}

		if(isScrambled) {
			return 1; // scrambled
		}

		thisAAs = new int[26];
		topAAs = new int[26];

		for(int i=0; i<thisStrip.length(); i++) {
			if(thisStrip.charAt(i) == 'E' || thisStrip.charAt(i) == 'D' || thisStrip.charAt(i) == 'Q' || thisStrip.charAt(i) == 'N') {
				continue;
			}
			thisAAs[thisStrip.charAt(i)-'A']++;
		}

		for(int i=0; i<thisTopStrip.length(); i++) {
			if(thisTopStrip.charAt(i) == 'E' || thisTopStrip.charAt(i) == 'D' || thisTopStrip.charAt(i) == 'Q' || thisTopStrip.charAt(i) == 'N') {
				continue;
			}
			topAAs[thisTopStrip.charAt(i)-'A']++;
		}

		boolean isDeamidated = true;
		for(int i=0; i<thisStrip.length(); i++) {
			if(thisAAs[i] != topAAs[i]) {
				isDeamidated = false;
			}
		}

		if(isDeamidated) {
			if(topPeptide.contains("+.98") || peptide.contains("+.98")) {
				return 2; // deamidation
			}
		}


		return 3; // others
	}
}

public class RerankedCategory {


	public static final int fractionIdx = 0;
	public static final int fileIdx = 1;
	public static final int scanIdx = 4;
	public static final int peptideIdx = 3;
	public static final int rankIdx = 18;
	public static final int isNoncanonicalIdx = 36;


	public static void main(String[] args) throws IOException {
		File peaksResult = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/6.CaliScanwoDeamiFix/PEAKS/PeptideAnnotationS4_5ppm_002_MSFrecal_addScanNum.tsv");
		File pxgResult = new File  ("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/6.CaliScanwoDeamiFix/pXg/PeptideAnnotationS4_5ppm_002_recal.scanNum.rep1.rank10.pXg.BA.fdr");

		Hashtable<String, String> a = spectrumToPeptide(peaksResult);
		Hashtable<String, PXGResult> b = spectrumToPeptidePXG(pxgResult);

		a.forEach((spec , peptide) -> {
			if(b.get(spec) != null) {
				b.get(spec).topPeptide = peptide;
			}
		});

		int[] cCategory = new int[4];
		int[] ncCategory = new int[4];
		String[] categories = {"Top-ranked", "Scrambled", "Deamidated", "Others"};

		b.forEach((spec , pxgRes) -> {
			int idx = pxgRes.getCategory();

			if(idx == 2) {
//				System.out.println(pxgRes.rank+"\t"+pxgRes.peptide+"\t"+pxgRes.topPeptide);
			}

			if(pxgRes.isCanonical.equalsIgnoreCase("TRUE")) {
				cCategory[idx]++;
			} else {
				ncCategory[idx]++;
			}

//			System.out.println(pxgRes.toString());
		});

		System.out.println("Category,Canonical PSMs,Noncanonical PSMs");
		for(int i=0; i<=3; i++) {
			System.out.println(categories[i]+","+cCategory[i]+","+ncCategory[i]);
		}

	}

	public static Hashtable<String, PXGResult> spectrumToPeptidePXG (File file) throws IOException {
		Hashtable<String, PXGResult> psmTable = new Hashtable<>();

		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;

		BR.readLine(); // skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String scanID = fields[fractionIdx]+"_"+fields[fileIdx] +"_" +fields[scanIdx];
			PXGResult res = new PXGResult();

			res.peptide = fields[peptideIdx];
			res.record = line;
			res.isCanonical = fields[isNoncanonicalIdx];
			res.rank = Integer.parseInt(fields[rankIdx]);

			psmTable.put(scanID, res);
		}

		BR.close();

		return psmTable;
	}

	/**
	 * For original PEAKS
	 *
	 * @param file
	 * @return
	 * @throws IOException
	 */
	public static Hashtable<String, String> spectrumToPeptide (File file) throws IOException {
		Hashtable<String, String> psmTable = new Hashtable<>();

		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;

		BR.readLine(); // skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String scanID = fields[fractionIdx]+"_"+fields[fileIdx] +"_" +fields[scanIdx];
			String strip = fields[peptideIdx].replaceAll("[+-.0123456789\\(\\)*]*", "");
			if(strip.length() >= 8 && strip.length() <= 15) {
				if(psmTable.get(scanID) == null) {
					psmTable.put(scanID, fields[peptideIdx]);
				}
			}

		}

		BR.close();

		return psmTable;
	}
}
