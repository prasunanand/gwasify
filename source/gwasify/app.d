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
  else if(generator == "details"){
    chr_details(file_name);
  }
  else{
    readfile_bgen(file_name);
  }
}

void annotation_generator(string file_name){
  auto pipe = pipeShell("gunzip -c " ~ file_name);
  File input = pipe.stdout;

  File anno_file =  File("snps.txt", "w");
  File error_file =  File("error_log.txt", "w");
  size_t index = 1;
  foreach(line; input.byLine){
    auto chr = to!string(line).split(" ")[0];
    // writeln(chr);
    //auto anno = chr[0].split(":");
    //if(anno.length == 4){
      anno_file.writeln(chr);
    //}
    //else{
    //error_file.writeln(index, " => ", anno);
    //}
  index++;
  }

}

void chr_details(string file_name){
  //auto pipe = pipeShell("gunzip -c " ~ file_name);
  File input = File(file_name);

  File chr_variants = File("variants.txt", "w");

  foreach(line; input.byLine){
    auto chr = to!string(line).split(",");
    //writeln(chr);
    auto chr_name = chr[0];
    auto positions = chr[1].split(":");
    chr_variants.writeln(line, "\t", positions[0], "\t", positions[1]);
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


void readfile_bgen(const string file_bgen) {
  writeln("entered ReadFile_bgen");
  
  File infile = File(file_bgen);

  writeln("ALL SET!");

  //The first four bytes
  uint bgen_snp_block_offset = infile.rawRead(new uint[1])[0];

  //The header block
  uint bgen_header_length = infile.rawRead(new uint[1])[0];
  writeln("bgen_header_length => ", bgen_header_length);
  assert(bgen_header_length <= bgen_snp_block_offset);
  bgen_snp_block_offset -= 4;
  uint bgen_nsnps = infile.rawRead(new uint[1])[0];
  writeln("No. of variant = > ", bgen_header_length);
  bgen_snp_block_offset -= 4;
  uint bgen_nsamples = infile.rawRead(new uint[1])[0];
  writeln("No. of samples = > ", bgen_nsamples);
  bgen_snp_block_offset-=4;
  char[] magic_chars = infile.rawRead(new char[4]);
  writeln(magic_chars);

  size_t ignore = bgen_header_length - 20; // check
  if(ignore != 0)
    infile.rawRead(new char[ignore]);
  bgen_snp_block_offset -= ignore;

  //BitArray bgen_flags = BitArray(32, cast(ulong*)infile.rawRead(new char[4]));
  uint bgen_flags = infile.rawRead(new uint[1])[0];
  bgen_snp_block_offset -= 4;

  uint CompressedSNPBlocks = (bgen_flags) & 3;
  writeln("CompressedSNPBlocks => ", CompressedSNPBlocks);
  uint layout = (bgen_flags & (15 << 2)) >> 2;
  writeln("layout =>", layout);
  uint sample_ids_presence = (bgen_flags & (1 << 31)) >> 31;
  writeln("sample_ids_presence => ", sample_ids_presence);
  uint LongIds = (bgen_flags) & 0x4;
  //writeln(LongIds);

  if (layout == 0) {
    writeln("This value is not supported");
    exit(0);
  }

  //infile.rawRead(new char[bgen_snp_block_offset]);
  writeln(bgen_snp_block_offset);

  // sample identifier block

  uint bgen_LSI =  infile.rawRead(new uint[1])[0];
  bgen_snp_block_offset -= 4;
  writeln("bgen_LSI => ", bgen_LSI);
  uint N = infile.rawRead(new uint[1])[0];
  bgen_snp_block_offset -= bgen_LSI;
  writeln("N => ", N);
  for(uint i = 0; i <N; i++){
    ushort bgen_LS1_length = infile.rawRead(new ushort[1])[0];
    //writeln("bgen_LS1_length => ", bgen_LS1_length);
    string bgen_LS1 = cast(string)infile.rawRead(new char[bgen_LS1_length]);
    writeln("bgen_LS1 =>", bgen_LS1);
  }

}
