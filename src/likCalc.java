//Avinash Ramu, Conrad Lab, WUSTL.
//Code to calculate the beta binomial genotype likelihoods.
//Update 12/04 - use the beta binomial for the RA, use binomial allele distribution for RR, AA.
//Usage - "java likCalc Counts_file 6.505 9.084 sampleID > op_f"

import java.io.*;
import java.lang.Math;
import java.util.*;
import flanagan.analysis.*;


class likCalc 
{	
	//Calculates the beta binomial likelihood
	public static double lik_bb(double alpha, double beta, int n, int k)
	{
	    while(k>500 || n>500) {
			k = k/2;
			n = n/2;
	    }
	   
	    double t1  = k + alpha;
		double t2  = n - k + beta;
		double num = Stat.beta(t1, t2);		
		double denom = Stat.beta(alpha, beta);
		double log_lik  = Math.log10(num) - Math.log10(denom);
		//System.out.println("t1 is " + t1 + "\tt2 is " + t2);
		//System.out.println("num is " + num + "\tdenom is " + denom);
		//System.out.println("log_lik is " + log_lik);
		return log_lik;
	}

	//Calculates the binomial likelihood
	public static double lik_bin(int n, int k)
	{
		double p = 0.5;
		double log_lik  = n * Math.log10(p);
		//System.out.println("log_lik is " + log_lik);
		return log_lik;
	}

	//Normalizes the three likelihoods
	public static double[] normalize_lik(double log_lik_RR, double log_lik_RA, double log_lik_AA)
	{
		double min_log_lik = Math.min(log_lik_RR, Math.min(log_lik_RA, log_lik_AA));
		double norm_log_lik_RR = (int)(log_lik_RR - min_log_lik + 0.499);
	  	double norm_log_lik_RA = (int)(log_lik_RA - min_log_lik + 0.499);
	  	double norm_log_lik_AA = (int)(log_lik_AA - min_log_lik + 0.499);
  			
		if(norm_log_lik_RR > 255)
			norm_log_lik_RR = 255;

		if(norm_log_lik_RA > 255)
			norm_log_lik_RA = 255;

		if(norm_log_lik_AA > 255)
			norm_log_lik_AA = 255;

		if(norm_log_lik_RR < 0)
			norm_log_lik_RR = 0;

		if(norm_log_lik_RA < 0)
			norm_log_lik_RA = 0;

		if(norm_log_lik_AA < 0)
			norm_log_lik_AA = 0;

		double return_liks[] = new double[3];
		return_liks[0] = norm_log_lik_RR;
		return_liks[1] = norm_log_lik_RA;
		return_liks[2] = norm_log_lik_AA;

		return return_liks;

	}

	// The file is of format "chr pos ref ref_ReadCount alt  alt_ReadCount"
	//ARGUMENTS - countfile, alpha, beta, sampleName
 	public static void main(String args[])
  	{
  		try{

  			double BASE_ERROR = 0.001;
  			String countFile = args[0];
  			Double RRalpha = Double.valueOf(args[1]); 
  			Double RRbeta = Double.valueOf(args[2]);
  			Double RAalpha = Double.valueOf(args[3]); 
  			Double RAbeta = Double.valueOf(args[4]);  			 
  			Double AAalpha = Double.valueOf(args[5]); 
  			Double AAbeta = Double.valueOf(args[6]); 
  			String sName = args[7];


  			double[] GT_lik = new double[10];
  			
			System.out.println("#Counts file " + args[0]);
  			System.out.println("#RRAlpha, RRBeta " + RRalpha + " " + RRbeta);
  			System.out.println("#RAAlpha, RABeta " + RAalpha + " " + RAbeta);
  			System.out.println("#AAAlpha, AABeta " + AAalpha + " " + AAbeta);

  			Hashtable<String, Integer> GT = new Hashtable<String, Integer>();
  			GT.put("AA", new Integer(0));
  			GT.put("AC", new Integer(1));
  			GT.put("AG", new Integer(2));
  			GT.put("AT", new Integer(3));
  			GT.put("CC", new Integer(4));
  			GT.put("CG", new Integer(5));
  			GT.put("CT", new Integer(6));
  			GT.put("GG", new Integer(7));
  			GT.put("GT", new Integer(8));
  			GT.put("TT", new Integer(9));

  			//System.out.println (GT.get("AA"));

  			FileInputStream fstream = new FileInputStream(countFile);
  			// Get the object of DataInputStream
 			DataInputStream in = new DataInputStream(fstream);
	 		BufferedReader br = new BufferedReader(new InputStreamReader(in));
	  		String strLine;
	  		//Read File Line By Line
	  		while ((strLine = br.readLine()) != null)   {
	  			// Print the content on the console
	  			//System.out.println (strLine);
	  			String[] tokens = strLine.split("\t");
	  			if(tokens[0].matches("#.*")) {
	  				//System.out.println (strLine);
	  				continue;
	  			}
	  			//System.out.println (strLine);
	  			//System.out.println ("Number of tokens: " + tokens.length);  
	  			int nTok = tokens.length;

	  			// set all GL's to 255 initial, just change RR, RA and AA liks
	  			for(int i = 0; i<10; i++)
  					GT_lik[i] = 255;


	  			String chr = tokens[0];
	  			Long pos = Long.valueOf(tokens[1]);

	  			char ref = tokens[2].charAt(0);
	  			int ref_RC = Integer.parseInt(tokens[3]);	  			  			
	  			
	  			String RR;
	  			String RA;
	  			String AA;
	
	  			//System.out.println (chr);
	  			//System.out.println (pos);
	  			char alt= 'N';
	  			int alt_RC;

	  			if(nTok > 4) {
	  				alt = tokens[4].charAt(0);
	  				alt_RC = Integer.parseInt(tokens[5]);
	  			}
	  			else { // this is the Ref Ref case, no ALT
	  				alt = 'N';
	  				alt_RC = 0;
	  			}


	  			RR = new StringBuilder().append(ref).append(ref).toString();
	  			RA = new StringBuilder().append((char)Math.min(ref,alt)).append((char)Math.max(ref,alt)).toString(); // right ordering of alleles
	  			AA = new StringBuilder().append(alt).append(alt).toString();

	  					
		


	  			//System.out.println (ref + "\t" + ref_RC);
	  			//System.out.println (alt + "\t" + alt_RC);

	  			
	  			//System.out.println (GT.get(RR));
 			
	  			//log_lik_bb = -0.3010; // binomial
	  			//System.out.println ("Beta bin log_lik " + log_lik_bb);

	  			
	  			/*System.out.println ("bin log_lik " + log_lik_bb);*/


	  			double norm_bb_liks[] = new double[3];
				double log_lik_bb_RA   =  -10 * lik_bb(RAalpha, RAbeta, ref_RC+alt_RC, alt_RC); // try sum of both liks
				double log_lik_bb_RR   =  -10 * lik_bb(RRalpha, RRbeta, ref_RC+alt_RC, alt_RC);
				double log_lik_bb_AA   =  -10 * lik_bb(AAalpha, AAbeta, ref_RC+alt_RC, alt_RC);
				System.out.println ("BB lik " + log_lik_bb_RR + "\t" + log_lik_bb_RA + "\t" + log_lik_bb_AA);
				norm_bb_liks =  normalize_lik(log_lik_bb_RR, log_lik_bb_RA, log_lik_bb_AA);
				log_lik_bb_RR = norm_bb_liks[0];
				log_lik_bb_RA = norm_bb_liks[1];
				log_lik_bb_AA = norm_bb_liks[2];
				System.out.println ("BB norm lik " + log_lik_bb_RR + "\t" + log_lik_bb_RA + "\t" + log_lik_bb_AA);




	  			/*double log_lik_RR = -1 * ( log_lik_bb + 
	  			                             ref_RC * Math.log10(2* (1-BASE_ERROR)) + 
	  			                             alt_RC * Math.log10(2* (BASE_ERROR)));// log for precision multiplication
	  			double log_lik_RA = -1 * log_lik_bb;
	  			double log_lik_AA = -1 * ( log_lik_bb + 
	  			                             alt_RC * Math.log10(2 * (1-BASE_ERROR)) + 
	  			                             ref_RC * Math.log10(2 * (BASE_ERROR)));
	  			System.out.println ("bin - log_lik " + log_lik_RR + " " +  log_lik_RA + " " + log_lik_AA);	  */
	  			double norm_binomial_liks[] = new double[3];
				double log_lik_binomial =  lik_bin(ref_RC+alt_RC, alt_RC); 
				double log_lik_RR_binomial = -10 * ( log_lik_binomial + 
	  			                             ref_RC * Math.log10(2* (1-BASE_ERROR)) + 
	  			                             alt_RC * Math.log10(2* (BASE_ERROR)));// log for precision multiplication
				double log_lik_RA_binomial = -10 * log_lik_binomial;
	  			double log_lik_AA_binomial = -10 * ( log_lik_binomial + 
	  			                             alt_RC * Math.log10(2 * (1-BASE_ERROR)) + 
	  			                             ref_RC * Math.log10(2 * (BASE_ERROR)));
				System.out.println ("binomial lik " + log_lik_RR_binomial + "\t" + log_lik_RA_binomial + "\t" + log_lik_AA_binomial);
	  			norm_binomial_liks =  normalize_lik(log_lik_RR_binomial, log_lik_RA_binomial, log_lik_AA_binomial);
				log_lik_RR_binomial = norm_binomial_liks[0];
				log_lik_RA_binomial = norm_binomial_liks[1];
				log_lik_AA_binomial = norm_binomial_liks[2];
				System.out.println ("binomial norm lik " + log_lik_RR_binomial + "\t" + log_lik_RA_binomial + "\t" + log_lik_AA_binomial);


				//System.out.println ("bin - log_lik " + norm_binomial_liks[0] + " " +  norm_binomial_liks[1] + " " + norm_binomial_liks[2]);	 
	  			

	  			//double norm_bb_liks[] = new double[3];
	  			//double log_lik_RR = /*log_lik_bb_RR; // + */ log_lik_RR_binomial; ///2; log for precision multiplication
	  		    //double log_lik_RA = /*log_lik_bb_RA; // + */ log_lik_RA_binomial; ///2;
	  			//double log_lik_AA = /*log_lik_bb_AA; // + */ log_lik_AA_binomial; ///2;
	  			
	  			double log_lik_RR = log_lik_bb_RR; 
	  		    double log_lik_RA = log_lik_bb_RA; 
	  			double log_lik_AA = log_lik_bb_AA; 
	  			//norm_bb_liks =  normalize_lik(log_lik_RR_bb, log_lik_RA_bb, log_lik_AA_bb);
	  			//System.out.println ("betbin - log_lik " + norm_bb_liks[0] + " " +  norm_bb_liks[1] + " " + norm_bb_liks[2]);	 


	  			//double norm_liks[] =  normalize_lik(log_lik_RR_binomial, log_lik_RA_bb, log_lik_AA_binomial);
	  			double norm_liks[] =  normalize_lik(log_lik_RR, log_lik_RA, log_lik_AA);
	  			double norm_log_lik_RR = norm_liks[0];
	  			double norm_log_lik_RA = norm_liks[1];
	  			double norm_log_lik_AA = norm_liks[2];
	  			System.out.println ("sum normalized log_lik " + norm_log_lik_RR + "\t" +  norm_log_lik_RA + "\t" + norm_log_lik_AA);

	  			//System.out.println ("bet + bin - log_lik " + log_lik_RR + " " +  log_lik_RA + " " + log_lik_AA);	  
	
	  			//double min_log_lik = Math.min(log_lik_RR, Math.min(log_lik_RA, log_lik_AA));
	  			//double max_log_lik = Math.max(log_lik_RR, Math.max(log_lik_RA, log_lik_AA));


	  			/*double norm_log_lik_RR = (log_lik_RR - min_log_lik)/(max_log_lik) * 255;
	  			double norm_log_lik_RA = (log_lik_RA - min_log_lik)/(max_log_lik) * 255;
	  			double norm_log_lik_AA = (log_lik_AA - min_log_lik)/(max_log_lik) * 255;
	  			

	  			double norm_log_lik_RR = log_lik_RR;
	  			double norm_log_lik_RA = log_lik_RA;
	  			double norm_log_lik_AA = log_lik_AA;*/

	  			// set all GL's to max, just change RR, RA and AA liks - exp cuatro
	  			//for(int i = 0; i<10; i++)
  					//GT_lik[i] = max_log_lik;
	  			/*double norm_log_lik_RR = (int)(log_lik_RR - min_log_lik + 0.499);
	  			double norm_log_lik_RA = (int)(log_lik_RA - min_log_lik + 0.499);
	  			double norm_log_lik_AA = (int)(log_lik_AA - min_log_lik + 0.499);


	  			
	  			if(norm_log_lik_RR > 255)
	  				norm_log_lik_RR = 255;

	  			if(norm_log_lik_RA > 255)
	  				norm_log_lik_RA = 255;

	  			if(norm_log_lik_AA > 255)
	  				norm_log_lik_AA = 255;

	  			if(norm_log_lik_RR < 0)
	  				norm_log_lik_RR = 0;

	  			if(norm_log_lik_RA < 0)
	  				norm_log_lik_RA = 0;

	  			if(norm_log_lik_AA < 0)
	  				norm_log_lik_AA = 0;*/
	  			
				//System.out.println ("log_lik " + norm_log_lik_RR + "\t" +  norm_log_lik_RA + "\t" + norm_log_lik_AA);

	  			
	  			if(nTok == 4) { // only reference allele.
	  				Enumeration e = GT.keys();
   
				    //iterate through Hashtable keys Enumeration
				    while(e.hasMoreElements()) {
				      //GT = new StringBuilder().append(e.nextElement().toString());
				   	  String GT1 = e.nextElement().toString();
				   	  int index = GT1.indexOf(ref);
				   	  if(index != -1) {
				      	//System.out.print(GT1 + " RA \t" + index + "\n");
				      	//GT_lik[GT.get(GT1)] = norm_log_lik_RA;// + log_lik_RA_binomial)/2; // GTs with a single Ref or more get RA lik value
				      	GT_lik[GT.get(GT1)] = 127; // 1211 check
				      }
				      else {
				      	//System.out.print(GT1 + " AA \t" + index + "\n");
				      	//GT_lik[GT.get(GT1)] = norm_log_lik_AA;
				      	GT_lik[GT.get(GT1)] = 255; // 1211 check
				      }

				  	}
				  	RR = new StringBuilder().append(ref).append(ref).toString(); // two copies or ref gets RR
	  				GT_lik[GT.get(RR)] = norm_log_lik_RR ;//+ log_lik_RR_binomial)/2;

	  				//System.out.print("\nREFREF " + RR + "\n");
	  				System.out.print (chr + "\t" + pos + "\t" + ref + "\t" + ref_RC + "\t" + alt + "\t" + alt_RC);
	  				for(int i = 0; i<10; i++)
  						System.out.print ("\t" + GT_lik[i]);
  					System.out.print("\t" +sName + "\n");
  
  					//System.exit(0);
  					continue;  					
  				}

				
  				GT_lik[GT.get(RR)] = norm_log_lik_RR;
	  			GT_lik[GT.get(RA)] = norm_log_lik_RA;
	  			GT_lik[GT.get(AA)] = norm_log_lik_AA;

	  			System.out.print (chr + "\t" + pos + "\t" + ref + "\t" + ref_RC + "\t" + alt + "\t" + alt_RC);
	  			for(int i = 0; i<10; i++)
  					System.out.print ("\t" + GT_lik[i]);
  				System.out.print("\t" + sName + "\n");



	  		}
	  		in.close();
    	}
    	catch (Exception e){
    		//Catch exception if any
  			System.err.println("Error: " + e.getMessage());
  		}
  	}
}