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

public class Landloop {

	public static void main(String[] args) {
		int N = 20;
		int K;
		long seed;
		String [] argv = new String[10];
		String directory = "lands/";
		argv[6] = "-o";
		argv[4] = "-s";
		argv[2] = "-k";
		argv[0] = "-n";
		argv[1] = Integer.toString(N);
		argv[8] = "-l";
		PrintStream out;
		try {
			out = new PrintStream(new File("files.txt"));
			//for(int n : N){
			for(int i = 0; i < 10; i++){
				for(int k = 0; k < N; k++){
					seed = (long)(Math.random()*50000);
					String file = directory +N+"-"+k+"-"+seed+".txt";
					argv[3] = Integer.toString(k);
					argv[5] = Long.toString(seed);
					argv[7] = file;
					argv[9] = directory + "landscape"+ N+"-"+k+"-"+seed+".txt";
					write(out,argv[9]);
					Configuration config = new Configuration(argv);
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
