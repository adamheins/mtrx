# mtrx
The mtrx project is a linear algebra library that contains common vector and matrix operations. Supported operations include matrix manipulation, matrix arithmetic, and linear system solving.

For a full list and description of all supported operations, view the [mtrx.h](https://github.com/adamheins/mtrx/blob/master/src/mtrx.h) and [vctr.h](https://github.com/adamheins/mtrx/blob/master/src/vctr.h) header files.

## Set Up
To use the library, start by cloning this repository. An example program can be found in the example/ folder. Run `make` to generate a binary called 'mtrx' that includes any c files in the example/ and src/ directories. Execute it by running `./mtrx`.

## Usage
If you only want the functionality of vectors in your program just write:
<pre>
#include &lt;sclr.h&gt;
#include &lt;vctr.h&gt;
</pre>

To gain matrices as well, write the following in addition to the above includes:
<pre>
#include &lt;mtrx.h&gt;
</pre>

## Testing
To run all of the tests, simply run `make test`. The source code for all the tests can be found in the tests/ folder.

This project uses clar for its unit testing framework. Source code and documentation for clar can be ofound [here](https://github.com/vmg/clar).

## Contributing
If you're interesting in improving the exising linear algebra algorithms or adding new ones, please feel free to submit a pull request!

## License
MIT license. For the full terms, please see the included LICENSE file.
