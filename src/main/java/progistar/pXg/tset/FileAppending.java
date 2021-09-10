package progistar.pXg.tset;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class FileAppending {

	/**
	 * Simply testing file appending...
	 * 
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		File file = new File("test.txt");
		
		BufferedWriter BW = new BufferedWriter(new FileWriter(file));
		BW.append("1");
		BW.close();
		
		BW = new BufferedWriter(new FileWriter(file, true));
		
		BW.append("2");
		BW.close();
	}
}
