# Setting up

To run this demo, you'll need to have
[Docker](//docs.docker.com/engine/install/) (minimum version 27) and
[Docker Compose](//docs.docker.com/compose/install/) (minimum version 2.17.0)
installed, as well as [`jq`](https://jqlang.github.io/jq/).

We also recommend a Node.js package manager (such as `npm` or `yarn`) to make
things a easier, but we will provide directions for doing everything without one
that you can follow if you prefer.

## Getting the data

First download the
[demo data](//s3.us-east-1.amazonaws.com/jbrowse.org/apollo/data.zip) used for
this guide and put it in a directory where you can also put a couple other
files. For example:

```sh
mkdir apollo-demo/
cd apollo-demo/
wget https://s3.us-east-1.amazonaws.com/jbrowse.org/apollo/data.zip
unzip data.zip
rm data.zip
```

You'll now have two directories called `data/` and `jbrowse_data/`.

## Running Apollo

Now in the terminal, run `docker compose up`. You should see a stream of logs
from the Docker containers. If you do, Apollo is now running! You can use
<kbd>Ctrl</kbd> + <kbd>C</kbd> in the terminal to stop Apollo when you are done.
