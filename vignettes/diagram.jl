# Draw a Graphviz diagram in a vignette. Needs the `dot` program (https://graphviz.org) on the
# PATH when rendering. Use in a cell with `#| echo: false`:
#
#     diagram("""
#     digraph { S -> I [label = "βSI"] }
#     """)
#
# The result is an SVG, so Quarto treats it like a plot in HTML, PDF and Markdown.

struct Diagram
    svg::String
end
Base.show(io::IO, ::MIME"image/svg+xml", d::Diagram) = print(io, d.svg)

const DIAGRAM_STYLE = """
    graph [rankdir = LR, fontname = "Helvetica", fontsize = 11, nodesep = 0.35, ranksep = 0.5,
           bgcolor = "transparent"];
    node  [shape = box, style = "rounded,filled", fillcolor = "#f2f2f2", color = "#555555",
           fontname = "Helvetica", fontsize = 11, margin = "0.12,0.06"];
    edge  [color = "#444444", fontname = "Helvetica", fontsize = 9, arrowsize = 0.7];
"""

"""
    diagram(dot)

Render Graphviz `dot` source to SVG. The house style (left to right, rounded boxes) is
inserted after the opening brace; set attributes in the source to override it.
"""
function diagram(dot::AbstractString)
    i = findfirst('{', dot)
    i === nothing && throw(ArgumentError("not a Graphviz graph"))
    src = dot[1:i] * "\n" * DIAGRAM_STYLE * dot[nextind(dot, i):end]
    svg = read(pipeline(IOBuffer(src), `dot -Tsvg`), String)
    # drop the XML prolog and doctype so the SVG can be embedded inline
    j = findfirst("<svg", svg)
    return Diagram(svg[first(j):end])
end
