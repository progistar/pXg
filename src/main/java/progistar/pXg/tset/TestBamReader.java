package progistar.pXg.tset;

import java.io.File;
import java.io.IOException;
import java.util.Map.Entry;
import java.util.Set;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class TestBamReader {

	public static void main(String[] args) {
        String bamFilePath = "/Users/gistar/eclipse-workspace/pXg/test/decoy_position_test/decoy_position.sam";

        String lineSeparator = System.getProperty("line.separator");
        try (SamReader samReader = SamReaderFactory.makeDefault().open(new File(bamFilePath))) {
        	
        	// get header
        	SAMFileHeader header = samReader.getFileHeader();
        	
        	Set<Entry<String, String>> attrs = header.getAttributes();
        	
        	attrs.forEach(s -> {
        		System.out.println(s.getKey());
        	});
        	
        	
            for (SAMRecord samRecord : samReader) {
                // Process each SAMRecord as needed
                System.out.println(samRecord.getSAMString().replace(lineSeparator, ""));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
