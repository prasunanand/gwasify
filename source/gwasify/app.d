module gwasify.app;

import std.getopt;
import std.stdio;
import std.string;

void main(string[] args){
  string file_name;

  getopt(args,
    "file"    , &file_name,
  );

  writeln(file_name);

  File pheno_covar_file = File(file_name);
  File pheno_file =  File("pheno.txt", "w");
  File covar_file =  File("covar.txt", "w");

  pheno_covar_file.readln();

  foreach(line; pheno_covar_file.byLine){
    auto items = line.split("\t");
    pheno_file.writeln(items[0]);

    covar_file.write(items[1]);
    foreach(covariate; items[2..$]){
      covar_file.write("\t");
      covar_file.write(covariate);
    }
    covar_file.write("\n");
  }
}

