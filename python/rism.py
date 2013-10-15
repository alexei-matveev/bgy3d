#!/usr/bin/python
# -*- mode: python; c-basic-offset: 4 -*- vim: set sw=4 tw=70 et sta ai:

import sys
import os
import shlex
from tempfile import mkdtemp
from subprocess import Popen
from contextlib import contextmanager


def sexp (x):
    "Convert array of numbers to a nested s-expression."
    try:
        n = len (x)
    except:
        return repr (x)
    return "(" + " ".join (sexp (y) for y in x) + ")"

#
# Request  is  an s-expression  in  text  form  encoding an  array  of
# numbers, reply is just a number again in text form:
#
def protocol (server):
    def f (x):
        return float (server (sexp (x)))
    return f


#
# The bulk of python code  operates with callable Objects of type Func
# for historical  reasons. So that also a  server-based solution needs
# to be cast  into that pattern. This is supposed to  be one such Func
# that happens to implement derivatives by finite differences:
#
@contextmanager
def Server (args):
    # FIXME: this  is convenience hack, one should  decide either args
    # is a list or a string:
    if type (args) == type (""):
        args = shlex.split (args)

    # Root process such as naive backup can still enter this directory
    # and read/write to the fifo pipes:
    tmp = mkdtemp()

    # Just one possible choice of names:
    inp = os.path.join (tmp, "%inp")
    out = os.path.join (tmp, "%out")

    os.mkfifo (inp)
    os.mkfifo (out)

    # args[]  is supposed  to include  the executable  and  flags, the
    # subcommand will be inserted. FIXME: this is ugly.
    proc = Popen ([args[0]] + ["server", inp, out] + args[1:])

    # This  function  takes  input  text  and returns  output  of  the
    # subprocess, both communicated via the fifos:
    def server (x):
        with open (inp, "w") as f:
            f.write (x)
        with open (out, "r") as f:
            y = f.read ()
        return y

    # The control returns to the "with" statement body:
    yield protocol (server)

    # When leaving  the "with" statement body execute  this.  Tell the
    # daemon process to terminate:
    with open (inp, "w") as f:
        f.write ("#f")          # convention

    proc.wait ()
    os.unlink (inp)
    os.unlink (out)
    os.rmdir (tmp)


def main (cmd):
    #
    # It does  not appear feasible to leak  the implementation details
    # of  the  server-based  solution   to  the  bulk  of  the  python
    # code. That code should still belive it works with a callable.
    #
    from pts.func import NumDiff
    from numpy import array

    x = array ([[-0.2929, 0.0, 0.0],
                [0.2929, 0.757, 0.0],
                [0.2929, -0.757, 0.0]])

    with Server (["/home/alexei/darcs/bgy3d/guile/runbgy.scm"]) as f:
        for _ in range (10):
            print f(x)
        g = NumDiff (f)
        print g.fprime(x)


if __name__ == "__main__":
    main (*sys.argv)
