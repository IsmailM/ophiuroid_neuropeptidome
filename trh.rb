require 'bio'

@opt = {
  input: 'trh_ordered1.fa',
  translated: 'trh.aa.fa'
}

refs = {}
ref = Bio::FlatFile.open(Bio::FastaFormat, @opt[:translated]).each do |entry|
  refs[entry.entry_id] = { clean_aa: entry.aaseq.gsub('-',''), aa: entry.aaseq }
end

Bio::FlatFile.open(Bio::FastaFormat, @opt[:input]).each do |entry|
  (1..6).each do |f|
    translated = entry.naseq.translate(f)
    next unless translated.include? refs[entry.entry_id][:clean_aa]
    translated.scan(refs[entry.entry_id][:clean_aa]) do |m|
      start_idx = Regexp.last_match.offset(0)[0] * 3
      end_idx   = (Regexp.last_match.offset(0)[1] * 3) - 1
      subset = entry.naseq[start_idx..end_idx]
      refs[entry.entry_id][:subset] = subset
      refs[entry.entry_id][:frame] = f
      refs[entry.entry_id][:transcript] = entry.naseq 
    end
  end
end

refs.each do |id, d|
  next if d[:subset].nil?
  gaps = []
  d[:aa].chars.each_with_index  { |c, i| gaps << i*3 if c == '-' }
  gaps.each { |g| d[:subset].insert(g, '---') }
  d[:aa] = ' ' + d[:aa].chars.join('  ')
  d[:subset] = d[:subset].upcase()
  puts ">#{id}"
  puts "#{d[:aa]}"
  puts "#{d[:subset]}"
end
