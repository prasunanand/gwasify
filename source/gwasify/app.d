module gwasify.app;

import std.getopt;
import std.stdio;
import std.string;
import core.stdc.stdlib : exit;
import std.conv;

void main(string[] args){
  string file_name;
  string pheno_cols, covar_cols;

  size_t mode = 1;
  getopt(args,
    "file"    , &file_name,
    "mode"  , &mode,
    "pheno_cols", &pheno_cols,
    "covar_cols", &covar_cols
  );

  auto pheno_arr = pheno_cols.split(",");
  auto covar_arr = covar_cols.split(",");

  writeln(file_name);

  File pheno_covar_file = File(file_name);
  File pheno_file =  File("pheno.txt", "w");
  File covar_file =  File("covar.txt", "w");

  auto header = pheno_covar_file.readln();
  if(mode == 0){
    auto header_names = header.strip.split("\t");

    foreach(index, item; header_names){
      writeln( index, " => ", item);
    }
    exit(0);
  }

  foreach(line; pheno_covar_file.byLine){
    auto items = line.split("\t");

    pheno_file.write(items[to!size_t(pheno_arr[0])]);
    foreach(col; pheno_arr[1..$]){
      pheno_file.write(items[to!size_t(col)]);
    }
    pheno_file.write("\n");


    covar_file.write(items[to!size_t(covar_arr[0])]);
    foreach(col; covar_arr[1..$]){
      covar_file.write("\t");
      covar_file.write(items[to!size_t(col)]);
    }
    covar_file.write("\n");
  }
}

