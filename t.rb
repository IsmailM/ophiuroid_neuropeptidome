require 'bio'
aligned_file = 'results1.fa'

def init(aligned_file)
  input_read = {}
  biofastafile = Bio::FlatFile.open(Bio::FastaFormat, aligned_file)
  biofastafile.each_entry do |entry|
    input_read[entry.entry_id] = {
      
    }
    puts entry unless guess_sequence_type(entry.seq) == :nucleotide
    # sp_out = run_signalp(entry.aaseq)
    # seq = format_output(id, entry.aaseq, sp_out)
    # input_read << seq
  end
  # produce_results(input_read)
end

def guess_sequence_type(seq)
  cleaned_sequence = seq.gsub(/[^A-Z]|[NX]/i, '')
  type = Bio::Sequence.new(cleaned_sequence).guess(0.9)
  (type == Bio::Sequence::NA) ? :nucleotide : :protein
end



init(aligned_file)