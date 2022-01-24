package progistar.pXg.utils;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class Logger {

	private static BufferedWriter BW = null;
	
	public static void append (String string) {
		try {
			BW.append(string);
		}catch (IOException ioe) {
			
		}
	}
	
	public static void newLine () {
		try {
			BW.newLine();
		}catch (IOException ioe) {
			
		}
	}
	
	public static void close () {
		try {
			BW.close();
		}catch (IOException ioe) {
			
		}
	}
	
	public static void create (String fileName) {
		try {
			BW = new BufferedWriter(new FileWriter(fileName));
		}catch (IOException ioe) {
			
		}
	}
}
