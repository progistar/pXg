package progistar.pXg.data;

import java.util.ArrayList;
import java.util.Hashtable;

public class Global {

	public static Hashtable<String, Integer> SEQUENCE_HASH = new Hashtable<>();
	public static ArrayList<String> SEQUENCE_ARRAYLIST = new ArrayList<>();

	public static Hashtable<String, String[]> EXON_LENGTH_HASH = new Hashtable<>();
	public static Hashtable<String, String[]> PERCENT_FULL_DIST_HASH = new Hashtable<>();
	public static Hashtable<String, String[]> PERCENT_EXON_DIST_HASH = new Hashtable<>();
	public static Hashtable<String, String[]> PERCENT_CDS_DIST_HASH = new Hashtable<>();
	public static Hashtable<String, String[]> FROM_START_DIST_HASH = new Hashtable<>();
	public static Hashtable<String, String[]> FROM_STOP_DIST_HASH = new Hashtable<>();
	//public static String[] SEQUENCE_ARRAY;

	public static void putExonLengths (String key, String[] exonLengths, boolean enforce) {
		if(EXON_LENGTH_HASH.get(key) == null || enforce) {
			EXON_LENGTH_HASH.put(key, exonLengths);
		}
	}
	public static void putPercentFullDist (String key, String[] pFullDist, boolean enforce) {
		if(PERCENT_FULL_DIST_HASH.get(key) == null || enforce) {
			PERCENT_FULL_DIST_HASH.put(key, pFullDist);
		}
	}
	public static void putPercentExonDist (String key, String[] pExonDist, boolean enforce) {
		if(PERCENT_EXON_DIST_HASH.get(key) == null || enforce) {
			PERCENT_EXON_DIST_HASH.put(key, pExonDist);
		}
	}
	public static void putPercentCDSDist (String key, String[] pCDSDist, boolean enforce) {
		if(PERCENT_CDS_DIST_HASH.get(key) == null || enforce) {
			PERCENT_CDS_DIST_HASH.put(key, pCDSDist);
		}
	}
	public static void putStartDist (String key, String[] startDist, boolean enforce) {
		if(FROM_START_DIST_HASH.get(key) == null || enforce) {
			FROM_START_DIST_HASH.put(key, startDist);
		}
	}
	public static void putStopDist (String key, String[] stopDist, boolean enforce) {
		if(FROM_STOP_DIST_HASH.get(key) == null || enforce) {
			FROM_STOP_DIST_HASH.put(key, stopDist);
		}
	}

	public static int getSequenceValue (String sequence) {
		Integer value = SEQUENCE_HASH.get(sequence);
		if(value == null) {
			value = SEQUENCE_HASH.size();
			SEQUENCE_HASH.put(sequence, value);
			SEQUENCE_ARRAYLIST.add(sequence);
		}

		return value;
	}

	/**
	 * @deprecated
	 */
	/*
	public static void updateSequenceArray () {
		long startTime = System.currentTimeMillis();
		System.out.println("Index sequence array...");

		SEQUENCE_ARRAY = new String[SEQUENCE_HASH.size()];

		SEQUENCE_HASH.forEach((sequence, value) -> {
			SEQUENCE_ARRAY[value] = sequence;
		});

		long endTime = System.currentTimeMillis();

		System.out.println("A total of "+SEQUENCE_ARRAY.length+" were indexed");
		System.out.println((endTime - startTime)/1000 +" sec");
	}
	*/

}
