# Draws etc/logo.png, in the style of the FinancialToolbox.jl logo: three discs in the Julia colours.
# Needs Luxor.jl in the active environment and the Georgia font: julia etc/logo.jl
using Luxor

const WIDTH, HEIGHT = 604, 552
const RING_RADIUS, DISC_RADIUS = 135, 117
const THIN_STROKE = 3
const BOLD_PATH_STROKE = 10
const SYMBOL_SCALE = 0.78 # keeps a margin between the symbols and the rings
const PATHS_SCALE = 0.9
const BLACK = "#000000"

# Ring and disc colours of the Julia logo.
const GREEN = ("#389826", "#60ad51")
const RED = ("#cb3c33", "#d5635c")
const PURPLE = ("#9558b2", "#aa79c1")

function draw_disc(centre, colours)
    setcolor(colours[1])
    circle(centre, RING_RADIUS, :fill)
    setcolor(colours[2])
    circle(centre, DISC_RADIUS, :fill)
end

function draw_line(points)
    setcolor(BLACK)
    poly(points, :stroke, close = false)
end

# Monte Carlo: paths leaving the same point, a few bold ones over the crowd.
function draw_paths()
    steps = -96:28:72
    thin_moves = ([-16, -6, -14, -4, -20, -6], [-6, -14, -4, -16, 2, -12], [-2, -12, 6, -14, -6, 0], [4, -10, -6, 10, -8, 2], [6, 8, -4, 14, -2, 8], [12, 2, 14, -6, 18, 12], [16, 10, 6, 14, -4, 26])
    bold_moves = ([-26, -22, 12, -30, 6, -22], [14, -24, 20, 10, -22, 14], [28, 14, 22, -8, 26, 2])
    setline(THIN_STROKE)
    foreach(move -> draw_line(Point.(steps, cumsum([0; move]))), thin_moves)
    setline(BOLD_PATH_STROKE)
    foreach(move -> draw_line(Point.(steps, cumsum([0; move]))), bold_moves)
end

# The last step: the density, shaded where a call pays, which is what the vibrato integrates.
function draw_density()
    base, height, width, strike = 66, 132, 38, 10
    bell(x) = base - height * exp(-x^2 / (2 * width^2))
    tail = strike:2:100
    setcolor(RED[1])
    poly([Point(strike, base); Point.(tail, bell.(tail)); Point(100, base)], :fill)
    setline(7)
    draw_line([Point(-104, base), Point(104, base)])
    setline(9)
    draw_line(Point.(-100:2:100, bell.(-100:2:100)))
    setline(5)
    draw_line([Point(strike, base), Point(strike, bell(strike))])
end

# The derivatives, of any order, that come out of the price.
function draw_partial()
    setcolor(BLACK)
    fontface("Georgia")
    fontsize(290)
    textoutlines("∂", Point(0, 0), :fill, halign = :center, valign = :middle)
end

function draw_symbol(draw, centre, symbol_scale)
    gsave()
    translate(centre)
    scale(symbol_scale)
    draw()
    grestore()
end

function draw_logo(path)
    Drawing(WIDTH, HEIGHT, path)
    background(1, 1, 1, 0)
    setlinecap("round")
    setlinejoin("round")
    top, left, right = Point(302, 145), Point(146, 415), Point(458, 415)
    for (centre, colours, draw, symbol_scale) in ((top, GREEN, draw_paths, PATHS_SCALE), (left, RED, draw_density, SYMBOL_SCALE), (right, PURPLE, draw_partial, SYMBOL_SCALE))
        draw_disc(centre, colours)
        draw_symbol(draw, centre, symbol_scale)
    end
    finish()
end

draw_logo(joinpath(@__DIR__, "logo.png"))
