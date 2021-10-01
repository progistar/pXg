package progistar.pXg.tset;

import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MDTag {
	private static final Pattern EACH_MD_REGEX = Pattern.compile("(([0-9]+)|([A-Z]+|\\^[A-Z]+))");
	
	public static void main(String[] args) throws IOException {
		String tag = "12G2G73^AGCAGAG2C2CT4";
		String tag2 = "MD:Z:12G12AG26C31";
		
		Matcher matcher = EACH_MD_REGEX.matcher(tag);
		
		while (matcher.find()) {
			System.out.println(matcher.group());
		}
	}
}
