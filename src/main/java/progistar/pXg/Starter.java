package progistar.pXg;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.constants.RunInfo;
import progistar.pXg.data.PIN;
import progistar.pXg.data.pXgRecord;
import progistar.pXg.data.parser.ParameterParser;
import progistar.pXg.data.parser.pXgParser;
import progistar.pXg.processor.Master;
import progistar.pXg.utils.Logger;

/**
 * Hello world!
 *
 */
public class Starter
{
    public static void main( String[] args )
    {

    	long startTime = System.currentTimeMillis();
    	// set this application id by start time.
    	// the start time is supposed to be an unique id.
    	// this unique id is used to name the temporary files.
    	Constants.UNIQUE_RUN_ID = startTime;

    	// fail to parse params
    	if(ParameterParser.parseParams(args) == -1) {
    		System.exit(1);
    	}

    	int samFileSize = Parameters.sequenceFilePaths.length;

    	for(int si=0; si<samFileSize; si++) {
    		Parameters.CURRENT_FILE_INDEX = si;
    		Master.ready(Parameters.genomicAnnotationFilePath, Parameters.sequenceFilePaths[Parameters.CURRENT_FILE_INDEX], Parameters.peptideFilePath);
        	Master.run();

        	long endTime = System.currentTimeMillis();

        	RunInfo.printProcessedChromosomes();
        	RunInfo.printFilterStat();
        	System.out.println("\tTotal elapsed time: "+((endTime-startTime)/1000) + " sec with "+RunInfo.totalProcessedPeptides +" peptides and " +RunInfo.totalProcessedReads+" reads");
        	System.out.println();

        	Logger.append("\tTotal elapsed time: "+((endTime-startTime)/1000) + " sec with "+RunInfo.totalProcessedPeptides +" peptides and " +RunInfo.totalProcessedReads+" reads");
        	Logger.newLine();
        	Logger.newLine();
    	}
    	Logger.close();

    	// merge pXg tmp outputs
    	ArrayList<pXgRecord>[] tmpOutputs = new ArrayList[samFileSize];
    	Parameters.isStringent = false; // allow all identified peptides
    	for(int i=0; i<samFileSize; i++) {
    		try {
    			File tmpOutput = new File(Parameters.tmpOutputFilePaths[i]);
				tmpOutputs[i] = pXgParser.parse(tmpOutput);
				// delete tmp file
				tmpOutput.delete();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
    	}

    	pXgParser.writeMergedResult(tmpOutputs, Parameters.outputFilePath);

		// parser to PIN
		PIN.parseOutput();

    }
}