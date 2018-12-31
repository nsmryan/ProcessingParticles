
int num_dots = 200;

float time = 0.0; //<>// //<>//

float velocity = 0.02;
float acceleration = 5.0;

float spread = 150;

Dot[] dots;

color orange = color(204, 104, 0);
color purple = color(153,50,204);

float dist = 60;
boolean randomizing = true;

Circle circle;
Line line;
Square square;
Triangle triangle;

int shape_index = 0;
int last_shape_index = 0;
Shape current_shape = Shape.LINE;
Shape last_shape = current_shape;

State current_state = State.RANDOM;

enum Shape {
  LINE, CIRCLE, SQUARE, TRIANGLE,
}

enum State {
  RANDOM, STRAIGHT, TIMED
}

void setup() {
  size(640, 480);
  
  circle = new Circle(new PVector(width/2, height/2), dist);
  
  line = new Line(new PVector((width/2) - dist, (height/2) - dist), new PVector((width/2) + dist, (height/2) - dist));
  
  square = new Square(new PVector(width/2, height/2), dist, dist);

  triangle = new Triangle(new PVector(width/2, height/2), dist, dist);

  dots = new Dot[num_dots];

  for (int i = 0; i < num_dots; i++) {
    dots[i] = new Dot(new PVector((width/2) + randomGaussian() * spread, (height/2) + randomGaussian() * spread), new PVector(1.0, 1.0));
  }
}

void draw() {  
  float ms = millis();
  float delta_time = (time - ms) / 1000.0;

  velocity = 0.05;
  
  float noise_scaler = 10.0;
  float noise_detail = 1.0;
  
  time = ms;

  background(0, 0, 0);
  fill(color(204, 104, 0));

  for (int i = 0; i < num_dots; i++) {

    // accelerate
    float t = i / (float)num_dots;
    PVector surface_point = new PVector(0, 0);
    if (current_shape == Shape.LINE) {
      // line
      surface_point = line.start.copy();
      PVector offset = line.end.copy();
      offset.mult(t);
      surface_point.mult(1 - t);
      surface_point.add(offset);
    } else if (current_shape == Shape.SQUARE) {
        // square
        Line edge;
        if (t < 0.25) {
          edge = square.low;
          t *= 4;
        } else if (t < 0.5) {
          edge = square.right;
          t = (t - 0.25) * 4;
        } else if (t < 0.75) { 
          edge = square.high;
          t = (t - 0.50) * 4;
        } else {
          edge = square.left;
          t = (t - 0.75) * 4;
        }

        surface_point = edge.start.copy();
        PVector offset = edge.end.copy();
        offset.mult(t);
        surface_point.mult(1 - t);
        surface_point.add(offset);
    } else if (current_shape == Shape.CIRCLE) {
      // circle
      surface_point = circle.center.copy();
      PVector offset = PVector.fromAngle(t * 2 * PI);
      offset.mult(circle.radius);
      surface_point.add(offset);
    } else {
      // Triangle
      surface_point = new PVector(0, 0);
      Line edge;

      if (t < 0.33) {
        edge = new Line(new PVector(triangle.center.x - triangle.width, triangle.center.y - triangle.height),
                        new PVector(triangle.center.x + triangle.width, triangle.center.y - triangle.height));
        t = t * 3.0;
      } else if (t < 0.66) {
        edge = new Line(new PVector(triangle.center.x + triangle.width, triangle.center.y - triangle.height),
                        new PVector(triangle.center.x,                  triangle.center.y + triangle.height));
        t = (t - 0.33) * 3;
      } else {
        edge = new Line(new PVector(triangle.center.x - triangle.width, triangle.center.y - triangle.height),
                        new PVector(triangle.center.x,                  triangle.center.y + triangle.height));
        t = (t - 0.66) * 3;
      }

      surface_point = edge.start.copy();
      PVector offset = edge.end.copy();
      offset.mult(t);
      surface_point.mult(1 - t);
      surface_point.add(offset);
    }

    noStroke();    
     
    // move
    dots[i].pos.lerp(surface_point, velocity);

    //fill(color(204, 104, 0));
    fill(lerpColor(orange, purple, 1.724 * noise(dots[i].pos.x, dots[i].pos.y)));

    // draw
    ellipse(dots[i].pos.x, dots[i].pos.y, 6, 6);
  }

  if (current_state == State.TIMED) {
    int slice_size = 1600;
    int slice_elements = 3;
    int slice_width = (int)(slice_size / slice_elements);
    int slice = (int)(((int)ms % slice_size) / slice_width);

    if (slice == 0) {
      current_shape = Shape.TRIANGLE;
    } else if (slice == 1) {
      current_shape = Shape.SQUARE;
    } else {
      current_shape = Shape.CIRCLE;
    }
  }
  
  
  if (current_shape != last_shape) {
    if (current_state == State.RANDOM) {
      for (int i = 0; i < num_dots; i++) {
        Dot tmp = dots[i];
        int swap_index = (int)random(i, num_dots);
        dots[i] = dots[swap_index];
        dots[swap_index] = tmp;
      }
    }
  }
  last_shape = current_shape;
}

void mousePressed() {
  last_shape = current_shape;
  if (mouseButton == LEFT) {
    switch (current_shape) {
      case LINE:
      current_shape = Shape.TRIANGLE;
      break;
    case TRIANGLE:
      current_shape = Shape.SQUARE;
      break;
    case SQUARE:
      current_shape = Shape.CIRCLE;
      break;
    case CIRCLE:
      current_shape = Shape.LINE;
      break;
    }
  } else if (mouseButton == RIGHT) {
    switch (current_state) {
      case RANDOM:
        current_state = State.STRAIGHT;
        break;
      case STRAIGHT:
        current_state = State.TIMED;
        break;
      case TIMED:
        current_state = State.RANDOM;
        break;
    }
  }
}

class Circle  {
  float radius;
  PVector center;
  
  Circle(PVector c, float r) {
    radius = r;
    center = c;
  }
  
  PVector towards_edge(PVector pos) {
    float dist = 0.0;
    PVector surface_point = center.copy();
   
    surface_point.sub(pos);    
    dist = surface_point.mag();    
    dist -= radius;    
    surface_point.normalize();    
    surface_point.mult(dist);
    
    return surface_point; 
  }
}

class Line {
  PVector start;
  PVector end;
  
  Line(PVector s, PVector e) {
    start = s;
    end = e;
  }
  
  PVector towards_edge(PVector pos) {
    PVector towards = new PVector(0,0);
    SquaredDistancePointToLineSegment(start, end, pos, towards);
    return towards;
  }
}

class Square {
  PVector center;

  float width;
  float height;

  Line low;
  Line high;
  Line right;
  Line left;
  
  Square(PVector c, float w, float h) {
    center = c;
    width = w;
    height = h;
    low = new Line(new PVector(center.x - width, center.y + width),
                   new PVector(center.x + width, center.y + width));

    high = new Line(new PVector(center.x + width, center.y + width),
                    new PVector(center.x + width, center.y - width));

    right = new Line(new PVector(center.x + width, center.y - width),
                     new PVector(center.x - width, center.y - width));

    left = new Line(new PVector(center.x - width, center.y - width),
                    new PVector(center.x - width, center.y + width));
  }
}

class Triangle {
  PVector center;
  float width;
  float height;

  Triangle(PVector c, float w, float h) {
    center = c;
    width = w;
    height = h;
  }
}

class Dot {
  PVector pos;
  PVector vel;

  Dot(PVector p, PVector v) {
    pos = p;
    vel = v;
  }
}

// calculate the squared distance of a point P to a line segment A-B
// and return the nearest line point S
float SquaredDistancePointToLineSegment(PVector A, PVector B, PVector P, PVector S)
{
  float vx = P.x-A.x,   vy = P.y-A.y;   // v = A->P
  float ux = B.x-A.x,   uy = B.y-A.y;   // u = A->B
  float det = vx*ux + vy*uy; 

  if (det <= 0)
  { // its outside the line segment near A
    S.set(A);
    return vx*vx + vy*vy;
  }
  float len = ux*ux + uy*uy;    // len = u^2
  if (det >= len)
  { // its outside the line segment near B
    S.set(B);
    return sq(B.x-P.x) + sq(B.y-P.y);  
  }
  // its near line segment between A and B
  float ex = ux / sqrt(len);    // e = u / |u^2|
  float ey = uy / sqrt(len);
  float f = ex * vx + ey * vy;  // f = e . v
  S.x = A.x + f * ex;           // S = A + f * e
  S.y = A.y + f * ey;

  return sq(ux*vy-uy*vx) / len;    // (u X v)^2 / len
}
