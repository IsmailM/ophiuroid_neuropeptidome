require 'csv'
require 'bio'
require 'tempfile'

# Extending the Bio:Sequence:AA class
class Bio::Sequence::AA
  # Returns an array of all possible Open Reading Frame. Assumes that the stop
  #   codon is characterised by a non-word character i.e. '*' (as used by the
  #   bioruby translation function). Utilises a lookahead regex that advances
  #   through the Bio::Sequence::AA object (or a string) by a single character
  #   at a time.
  def findorfs(minsize = 10)
    scan(/(?=(M\w{#{minsize},}))./).flatten
  end

  def stop2stop(minsize = 10)
    scan(/\w{#{minsize},}/).flatten
  end
end

def parse_csv(blast_content, tab_format)
  parsed_csv = []
  lines = CSV.parse(blast_content, col_sep: ' ', skip_lines: /^#/,
                                   headers: tab_format)
  lines.each do |line|
    parsed_csv << line.to_hash
  end
  parsed_csv
end

def split_based_on_query(parsed_csv)
  flattened_csv = {}
  parsed_csv.each do |row|
    flattened_csv[row['sseqid']] ||= []
    flattened_csv[row['sseqid']] << row
  end
  flattened_csv
end

def analyse(flattened_csv)
  results = {}
  flattened_csv.each do |type, array|
    type_array = []
    array.each do |query|
      taxon     = format_taxon_name(query['taxon'])
      id        = '>' + taxon
      seq       = Bio::Sequence::NA.new(query['contig'])
      proteins  = translate(seq, query['frames'])
      query     = query['qseq'].gsub('-', '').gsub('*', '\*')
      full_orf  = proteins.scan(/\w*#{query}\w*/)
      idx       = ordered_number(taxon)
      type_array << { idx: idx, id: id, full_orf: full_orf }
    end
    results[type] = type_array.sort_by { |e| e[:idx] }
  end
  results
end

def ordered_number(taxon)
  order = %w(Am_cipu Am_cten Am_cons_1 Am_cons_2 Op_pilo Op_filo Op_wurd
             Mi_grac Am_squa Op_resi Op_savi Op_abys Op_angu Op_caes Op_fasc
             Op_scha Op_reti Op_lame Op_bisc Op_impr Op_bisp Op_vall Op_brev
             Ba_hero Op_tube Op_appr Op_macu Op_vivi Op_aust Op_cyli Op_wend
             Op_plic Op_obst Op_fune Op_perf Cl_cana Op_exim_1 Op_exim_2
             Op_liod Op_prol As_tubi As_bidw Op_oedi Go_pust As_love Op_iner
             Gl_sp_no Op_sp_no Am_laud Am_spat Op_john Op_lyma)
  # Ap_wils Lu_sene Le_tenu
  order.index(taxon)
end

def format_taxon_name(taxon)
  taxon.gsub!('.', '_')
  taxon.gsub!('_sp_cf_', '_')
  if taxon =~ /^\w+\d/
    taxon.scan(/^(..)\w+_(.. ?..).+(\d)/) do |a, e, d|
      return "#{a}_#{e.gsub(' ', '-')}_#{d}"
    end
  else
    taxon.scan(/^(..)\w+_(.._?..)/) { |a, e| return "#{a}_#{e.gsub(' ', '-')}" }
  end
end

def translate(sequence, _frame)
  seq = ''
  (1..6).each { |f| seq += "**||**#{sequence.translate(f)}" }
  seq
end

def stop2stop(seq)
  seq.stop2stop
end

def write_to_output(results)
  results.each do |type, array|
    uniq_array = ensure_uniq_fasta(array)
    write_results(type, uniq_array)
  end
end

def ensure_uniq_fasta(array)
  seqs = []
  array.delete_if do |e|
    exists = (seqs.include?(e[:full_orf]) || e[:full_orf].nil?)
    seqs << e[:full_orf] unless exists
    exists
  end
  array
end

def write_results(type, array)
  file = File.join(File.expand_path('original_fasta'),
                   "#{type.downcase.gsub('/', '_')}.fa")
  File.open(file, 'w') do |f|
    array.each do |e|
      f.puts e[:id]
      f.puts e[:full_orf]
    end
  end
end

@sp_path = '/Volumes/Data/data/programs/signalp-4.1/signalp'

blast_tab = ARGV[0]
blast_content = File.read(blast_tab)
tab_format = 'taxon qseqid sseqid qlen slen frames length nident positive' \
             ' mismatch sstart send score bitscore evalue qseq contig'

parsed_csv = parse_csv(blast_content, tab_format.split)
flattened_csv = split_based_on_query(parsed_csv)
results = analyse(flattened_csv)
write_to_output(results)
