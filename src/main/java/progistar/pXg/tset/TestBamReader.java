package progistar.pXg.tset;

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class TestBamReader {

	public static void main(String[] args) {
        String bamFilePath = "/Users/gistar/eclipse-workspace/pXg/test/decoy_position_test/decoy_position_.sam";

        try (SamReader samReader = SamReaderFactory.makeDefault().open(new File(bamFilePath))) {
        	
        	// get header
        	SAMFileHeader header = samReader.getFileHeader();
        	if (header.getTextHeader() != null) {
        		for (String line : header.getTextHeader().split("\n")) {
        			System.out.println(line);
        		}
        		
        	}
        	
        	// get comment
        	for (String comment : header.getComments()) {
        		System.out.println(comment);
        	}
        	// get records
            for (SAMRecord samRecord : samReader) {
                // Process each SAMRecord as needed
                System.out.print(samRecord.getSAMString());
                break;
            }
            
            for (SAMRecord samRecord : samReader) {
                // Process each SAMRecord as needed
                System.out.print(samRecord.getSAMString());
                break;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}