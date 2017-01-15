/**
 * Copyright (C) 2017, Sonia Singhal
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
package directional;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

/**
 * @author Sonia
 *
 */
public class SimLoop {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int N = 20;
		double mu = 0.05;
		int K;
		long seed;
		// double cuts[] = new double [] {0.1, 10, 0.71, 100};
		int gens = 25;
		String [] argv = new String[20];
		String directory = "outs-2Jan17-2/";
		argv[0] = "-n";
		argv[1] = Integer.toString(N);
		argv[2] = "-k";
		argv[4] = "-s";
		argv[6] = "-o";
		argv[8] = "-l";
		argv[10] = "-c";
		argv[11] = "{0.1 10 0.71 25}";
		// System.out.println(argv[11]);
		argv[12] = "-r";
		argv[13] = Double.toString(mu);
		argv[14] = "-g";
		argv[15] = Integer.toString(gens);
		argv[16] = "-f";
		argv[17] = "resources/population.txt";
		argv[18] = "-m";
		argv[19] = "multi_random";
		PrintStream out;
		try {
			out = new PrintStream(new File("files-2Jan17-3.txt"));
			//for(int n : N){
			for(int i = 0; i < 7; i++){
				for(int k = 0; k < N; k++){
					seed = (long)(Math.random()*50000);
					String file = directory +N+"-"+k+"-"+seed+".txt";
					argv[3] = Integer.toString(k);
					argv[5] = Long.toString(seed);
					argv[7] = file;
					argv[9] = directory + N+"-"+k+"-"+seed+".txt";
					write(out,argv[9]);
					// Configuration config = new Configuration(argv);
					// Simulation sim = new Simulation(argv);
					//sim.main(argv);
					Simulation.main(argv);
					System.out.println("done "+file);
				}

			}
			out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	static void write(PrintStream out,String line){
		out.println(line);
		return;
	}
}
