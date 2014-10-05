function usage()
io.stderr:write(string.format("Usage: %s <GFF file>\n" , arg[0]))
os.exit(1)
end

if #arg < 1 then
  usage()
  os.exit(1)
end

elem_out_v = gt.custom_visitor_new()
function elem_out_v:visit_feature(tn)
  local strand = "+"
  if tn:get_strand() == '-' then
    strand = '-'
  end
  for fn in tn:children() do
    if fn:get_type() == 'LTR_retrotransposon' then
      io.write("chr" .. fn:get_seqid() .. "\t"
                           .. fn:get_range():get_start()
                           .. "\t" .. fn:get_range():get_end()
                           .. "\telement"
                           .. "\t1\t" .. strand)
      local ltrpairs = {}
      for n in fn:children() do
        if n:get_type() == 'long_terminal_repeat' then
          table.insert(ltrpairs, {len=n:get_range():length(),
             start=(n:get_range():get_start()-fn:get_range():get_start())})
        end
      end
      io.write("\t0\t0\t100,100,100\t")
      io.write(#ltrpairs)
      io.write("\t")
      -- block lengths
      for i,item in ipairs(ltrpairs) do
        io.write(item.len-1)
        if i < #ltrpairs then
          io.write(',')
        end
      end
      io.write("\t")
      -- block starts
      for i,item in ipairs(ltrpairs) do
        io.write(item.start)
        if i < #ltrpairs then
          io.write(',')
        end
      end
      io.write("\n")
    end
  end
  return 0
end

dom_out_v = gt.custom_visitor_new()
function dom_out_v:visit_feature(fn)
  local strand = "+"
  if fn:get_strand() == '-' then
    strand = '-'
  end
  for n in fn:children() do
    if n:get_type() == 'protein_match' then
      print("chr" .. n:get_seqid() .. "\t" .. n:get_range():get_start()
                  .. "\t" .. n:get_range():get_end()
                 .. "\t" .. n:get_attribute("name")
                 .. "\t1\t" .. strand)
    end
    if n:get_type() == 'RR_tract' then
      print("chr" .. n:get_seqid() .. "\t" .. n:get_range():get_start()
                  .. "\t" .. n:get_range():get_end()
                 .. "\tPPT"
                 .. "\t1\t" .. strand)
    end
    if n:get_type() == 'primer_binding_site' then
      print("chr" .. n:get_seqid() .. "\t" .. n:get_range():get_start()
                  .. "\t" .. n:get_range():get_end()
                 .. "\tPBS"
                 .. "\t1\t" .. strand)
    end
  end
  return 0
end

vstream = gt.custom_stream_new_unsorted()
function vstream:next_tree()
  local node = self.instream:next_tree()
  if node then
     node:accept(self.visitor)
  end
  return node
end

vstream.instream = gt.gff3_in_stream_new_sorted(arg[1])
vstream.visitor = elem_out_v
print("track name=LTR description=\"LTRharvest candidates\" visibility=3")
local gn = vstream:next_tree()
while (gn) do
  gn = vstream:next_tree()
end

vstream.instream = gt.gff3_in_stream_new_sorted(arg[1])
vstream.visitor = dom_out_v
print("track name=protdoms description=\"LTRdigest features\" visibility=3")
local gn = vstream:next_tree()
while (gn) do
  gn = vstream:next_tree()
end

