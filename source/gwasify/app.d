module gwasify.app;

import core.stdc.stdlib : exit;

import std.conv;
import std.getopt;
import std.process;
import std.stdio;
import std.string;
import std.zlib;

void main(string[] args){
  string file_name;
  size_t mode = 1;
  string pheno_cols, covar_cols;
  string generator;

  getopt(args,
    "file"      , &file_name,
    "mode"      , &mode,
    "pheno_cols", &pheno_cols,
    "covar_cols", &covar_cols,
    "gen"       , &generator
  );

  if(generator == "covar"){
    pheno_covar_generator(file_name, mode, pheno_cols, covar_cols);
  }
  else if(generator == "snps"){
    annotation_generator(file_name);
  }
}

void annotation_generator(string file_name){
  auto pipe = pipeShell("gunzip -c " ~ file_name);
  File input = pipe.stdout;

  File anno_file =  File("snps.txt", "w");
  foreach(line; input.byLine){
    auto chr = to!string(line).split(",")[0];
    auto anno = chr.split(":");
    anno_file.writeln(anno[0], "\t", anno[1], "\t", anno[2], anno[3]);
  }

}

void pheno_covar_generator(string file_name, size_t mode, string pheno_cols, string covar_cols){
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
