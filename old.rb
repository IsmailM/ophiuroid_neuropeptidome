require 'bio'
require 'csv'
require 'haml'
require 'tempfile'

def init(aligned_file)
  input_read = []
  biofastafile = Bio::FlatFile.open(Bio::FastaFormat, aligned_file)
  biofastafile.each_entry do |entry|
    id = format_id(entry.entry_id)
    sp_out = run_signalp(entry.aaseq)
    seq = format_output(id, entry.aaseq, sp_out)
    input_read << seq
  end
  produce_results(input_read)
end

def format_id(id)
  id.scan(/^(..)\w+_(....)/) { |a, e| return "#{a}_#{e}" }
end

def format_output(id, seq, sp_out, seq_line_lengh = 95) #58
  a = []
  seq.scan(/.{1,#{seq_line_lengh}}/) do |s|
    index = Regexp.last_match.offset(0).first
    if sp_out['Sp?'] == 'Y'
      sp_start = convert_signalp_pos(seq, 1)
      sp_start_index = (sp_start / seq_line_lengh).ceil
      sp_end = convert_signalp_pos(seq, sp_out['Cmax-pos'])
      sp_end += 22
      sp_end_index = (sp_end / seq_line_lengh).ceil
      if index == sp_start_index
        s = s.insert(sp_start, '<span class="signalp">')
      end
      if index == sp_end_index
        s = s.insert(sp_end, '</span>')
        s.gsub!(/K-*K|K-*R|R-*R/) do |m|
          if Regexp.last_match.offset(0)[0] > sp_end
            "<span class=\"motif\">#{m}</span>"
          else
            m
          end
        end
      end
      if index != sp_start_index && index != sp_end_index
        s.gsub!(/K-*K|K-*R|R-*R/, '<span class="motif">\0</span>')
      end
    else
      s.gsub!(/K-*K|K-*R|R-*R/, '<span class="motif">\0</span>')
    end
    l =  "<b>#{id}</b>"
    l << '&nbsp;&nbsp;&nbsp;&nbsp;' if id.length == 7
    l << '&nbsp;&nbsp;&nbsp;' if id.length == 8
    l << '&nbsp;&nbsp;' if id.length == 9
    l << "#{s}"
    a << l
  end
  a
end

def produce_results(input_read)
  results = []
  input_read.each_with_index do |_, i|
    a = []
    input_read.map { |r| a << r[i] unless r[i].nil? }
    results << a unless a.empty?
  end
  results
end

def convert_signalp_pos(seq, sp_pos)
  pos = 0
  s = seq.chars
  s.each_with_index do |c, i|
    pos += 1 unless c == '-'
    return i if pos == sp_pos.to_i
  end
end

def colour_signalp(sp)
  seq = ''
  pos = 0
  s.each_char do |c|
    pos += 1 unless c == '-'
    seq << '<span class="signalp">' if pos == 1
    seq << c
    seq << '</span>' if pos == sp
  end
  seq
end

def signalp(input_file)
  sp_headers = %w(name Cmax Cmax-pos Ymax Ymax-pos Smax Smax-pos Smean D Sp?
                  Dmaxcut Networks-used)
  cmd = "#{@signalp} -t euk -f short -U 0.34 -u 0.34 #{input_file}"
  output = `#{cmd}`.gsub(/ +/, ',')
  sp_out = {}
  lines = CSV.parse(output, col_sep: ',', skip_lines: /^#/, headers: sp_headers)
  lines.each { |line| sp_out[line[0]] = line.to_hash }
  sp_out
end

def write_html_results(output, array)
  haml_doc = <<EOT
!!!
%html
  %head
    %title Results
    %meta{"http-equiv" => "Content-Type", :content => "text/html; charset=utf-8"}
    :css
      .id {font-weight: bold;}
      .signalp {background-color:#00f; color: #fff; font-weight: bold;}
      .motif {background-color:#29FD2E; color:#000; font-weight: bold;}
      p {word-wrap: break-word; font-weight: bold; font-size: 11pt; font-family:Courier New, Courier, Mono;}
  %body
    - array.each do |section|
      %p
        %span= section.join('</span><br><span>')
EOT
  engine = Haml::Engine.new(haml_doc)
  File.open(output, 'w') do |f|
    f.puts engine.render(Object.new, array: array)
  end
end

def run_signalp(contents)
  sp_headers = %w(name Cmax Cmax-pos Ymax Ymax-pos Smax Smax-pos Smean D Sp?
                  Dmaxcut Networks-used)
  file = Tempfile.new('signalp')
  file.write(">seq\n#{contents}")
  file.close
  cmd = "#{@signalp} -t euk -f short -U 0.34 -u 0.34 #{file.path} |" \
        " sed -n '3 p'"
  sp_out = `#{cmd}`
  Hash[sp_headers.zip(sp_out.split)]
ensure
  file.unlink
end

def run_muscle(input_file)
  output_file = "#{input_file}.muscle"
  cmd = "muscle -in #{input_file} -out #{output_file}"
  `#{cmd}`
  output_file
end

def run_mafft(input_file)
  output_file = "#{input_file}.mafft"
  cmd = "mafft --maxiterate 1000 --thread 8 #{input_file} > #{output_file}"
  `#{cmd}`
  output_file
end

@signalp = '/Volumes/Data/programs/signalp-4.1/signalp'
fasta_file = ARGV[0]
# aligned_file = run_muscle(fasta_file)
aligned_file = run_mafft(fasta_file)
output_file = "#{aligned_file}.doc"
input_read = init(aligned_file)
write_html_results(output_file, input_read)
