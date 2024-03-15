package progistar.pXg.data;

import java.util.Hashtable;

public class Global {

	public static Hashtable<String, Integer> SEQUENCE_HASH = new Hashtable<String, Integer>();
	public static String[] SEQUENCE_ARRAY;
	
	
	public static int getSequenceValue (String sequence) {
		Integer value = SEQUENCE_HASH.get(sequence);
		if(value == null) {
			value = SEQUENCE_HASH.size();
			SEQUENCE_HASH.put(sequence, value);
		}
		
		return value;
	}
	
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
	
}
