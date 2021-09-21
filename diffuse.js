var grid; //initial grid
var next; //updated grid
var vA = 1; //diffusion rate of the voltage

function setup() {
  createCanvas(200, 200);
  pixelDensity(1);
  grid = []; //declare grid as an array
  next = []; //declare next as an array

  //declare each element in the arrays as an array, creating a 2D array
  for (var x = 0; x < width; x++) { 
    grid[x] = [];
    next[x] = [];

    //assign the object to each array in the 2D array
    for (var y = 0; y < height; y++) {
      grid[x][y] = { v: 0 };
      next[x][y] = { v: 0 };
    }
  }

  //seed the voltage
  for (var j = 0; j < height; j++) {
    grid[1][j].v = 1;
  }
}

function draw() {
  background(0);

  //loop to access each cell
  for (var x = 1; x < width - 1; x++) {
    for (var y = 1; y < height - 1; y++) {
      //get the values from the object in each cell
      var v = grid[x][y].v;
      
      //update the values in the object using the algorithm and assign it to other array
      next[x][y].v = v + vA * laplaceV(x, y);
      //constrain values
      next[x][y].v = constrain(next[x][y].v, 0, 1);
    }
  }

  loadPixels();
  for (var x = 0; x < width; x++) {
    for (var y = 0; y < height; y++) {
      var pix = (x + y * width) * 4;
      var v = next[x][y].v;
      var c = floor(v * 255);
      c = constrain(c, 0, 255);
      pixels[pix + 0] = c;
      pixels[pix + 1] = c;
      pixels[pix + 2] = c;
      pixels[pix + 3] = 255;
    }
  }
  updatePixels();
  swap();
}

function laplaceV(x, y) {
  var sumV = 0;
  sumV += grid[x][y].v * -1;
  sumV += grid[x - 1][y].v * 0.67;
  sumV += grid[x][y].v * 0.33;
  return sumV;
}

function swap() {
  var temp = grid;
  grid = next;
  next = temp;
}
