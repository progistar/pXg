package progistar.pXg.utils;

import java.util.Hashtable;

public class IndexConvertor {

	private static Hashtable<String, Integer> strToIndex = new Hashtable<String, Integer>();
	private static Hashtable<Integer, String> indexToStr = new Hashtable<Integer, String>();
	
	/**
	 * Supported index: chr1, chr2, ..., chr22, chrX, chrY, chrMT. <br>
	 * chr1 to 0, chr2 to 1, ..., chrX to 22, chrY to 23, chrMT to 24. <br>
	 * If you give an unsupported index, it will return -1.
	 * 
	 * @param chr
	 * @return
	 */
	public static int chrToIndex (String chr) {
		int indexValue = -1;
		try {
			chr = chr.toLowerCase();
			indexValue = strToIndex.get(chr);
		} catch (Exception e) {
			
		}
		
		return indexValue;
	}
	/**
	 * 
	 * the index for that chr is automatically assigned by auto-increment key. <br>
	 * 
	 * @param chr
	 */
	public static void putChrIndexer (String chr) {
		chr = chr.toLowerCase();
		if(strToIndex.get(chr) != null) return;
		
		strToIndex.put(chr, strToIndex.size());
		indexToStr.put(indexToStr.size(), chr);
		
		if(strToIndex.size() != indexToStr.size()) {
			System.out.println(strToIndex.size() +":" +indexToStr.size());
			System.out.println("chr indexer has an error..!");
			System.out.println("::DO NOT MATCH INDEX!");
		}
	}
	
	public static String indexToChr (int index) {
		return indexToStr.get(index);
	}
}
