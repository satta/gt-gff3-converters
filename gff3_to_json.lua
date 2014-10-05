function usage()
io.stderr:write(string.format("Usage: %s <GFF file>\n" , arg[0]))
os.exit(1)
end

if #arg < 1 then
  usage()
  os.exit(1)
end

function json_encode(s)
  s = string.gsub(s,'\\','\\\\')
  s = string.gsub(s,'"','\\"')
  s = string.gsub(s,"'","\\'")
  s = string.gsub(s,'\n','\\n')
  s = string.gsub(s,'\t','\\t')
  return s
end

cv = gt.custom_visitor_new()
cv.first = true
function cv:print_feature_node(fn)
  if not cv.first then
    io.write(",")
  else
    cv.first = false
  end
  local score = "."
  if fn:get_score() then
    score = fn:get_score()
  end
  local attribs = {}
  io.write("\n")
  io.write("  {\n")
  io.write("    \"node\" : \"feature\",\n")
  io.write("    \"type\" : \"" .. json_encode(fn:get_type()) .. "\",\n")
  io.write("    \"source\" : \"" .. json_encode(fn:get_source()) .. "\",\n")
  io.write("    \"start\" : " .. fn:get_range():get_start() .. ",\n")
  io.write("    \"end\" : " .. fn:get_range():get_end() .. ",\n")
  io.write("    \"score\" : \"" .. score .. "\",\n")
  io.write("    \"strand\" : \"" .. fn:get_strand() .. "\",\n")
  io.write("    \"phase\" : \"" .. fn:get_phase() .. "\",\n")
  io.write("    \"attributes\" : {\n")
  for k,v in fn:attribute_pairs() do
    table.insert(attribs, "      \"" .. json_encode(k)
                              .. "\" : \"" .. json_encode(v) .. "\"")
  end
  io.write(table.concat(attribs, ",\n"))
  io.write("\n    }\n")
  io.write("  }")
  return 0
end
function cv:visit_feature(fn)
  for n in fn:children() do
    self:print_feature_node(n)
  end
  return 0
end
function cv:visit_region(n)
  if not cv.first then
    io.write(",")
  else
    cv.first = false
  end
  io.write("\n")
  io.write("  {\n")
  io.write("    \"node\" : \"region\",\n")
  io.write("    \"seqid\" : \"" .. json_encode(n:get_seqid()) .. "\",\n")
  io.write("    \"start\" : " .. n:get_range():get_start() .. ",\n")
  io.write("    \"end\" : " .. n:get_range():get_end() .. "\n")
  io.write("  }")
  return 0
end
function cv:visit_comment(n)
  if not cv.first then
    io.write(",")
  else
    cv.first = false
  end
  io.write("\n")
  io.write("  {\n")
  io.write("    \"node\" : \"comment\",\n")
  io.write("    \"comment\" : \"" .. json_encode(n:get_comment()) .. "\",\n")
  io.write("  }")
  return 0
end
function cv:visit_meta(n)
  if not cv.first then
    io.write(",")
  else
    cv.first = false
  end
  io.write("\n")
  io.write("  {\n")
  io.write("    \"node\" : \"meta\",\n")
  io.write("    \"directive\" : \"" .. json_encode(n:get_directive()) .. "\",\n")
  io.write("    \"data\" : \"" .. json_encode(n:get_data()) .. "\",\n")
  io.write("  }")
  return 0
end
function cv:visit_sequence(n)
  if not cv.first then
    io.write(",")
  else
    cv.first = false
  end
  io.write("\n")
  io.write("  {\n")
  io.write("    \"node\" : \"sequence\",\n")
  io.write("    \"header\" : \"" .. json_encode(n:get_description()) .. "\",\n")
  io.write("    \"sequence\" : \"" .. json_encode(n:get_sequence()) .. "\",\n")
  io.write("  }")
  return 0
end

json_out_stream = gt.custom_stream_new_unsorted()
json_out_stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
json_out_stream.visitor = cv
function json_out_stream:next_tree()
  local node = self.instream:next_tree()
  if node then
     node:accept(self.visitor)
  end
  return node
end

io.write("[")
out_stream = json_out_stream
local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end
io.write("\n]\n")
