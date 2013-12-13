#!/usr/bin/python
# -*- mode: python; c-basic-offset: 4 -*- vim: set sw=4 tw=70 et sta ai:
from __future__ import with_statement

import sys
import os
import shlex
from tempfile import mkdtemp
from subprocess import Popen
from contextlib import contextmanager
from pyparsing import Regex, Suppress, Forward, ZeroOrMore, Group
from numpy import array
from pts.func import Func
from pts.units import kcal


def to_sexp (x):
    "Convert array of numbers to a nested s-expression."
    try:
        n = len (x)
    except:
        return repr (x)
    return "(" + " ".join (to_sexp (y) for y in x) + ")"


def make_sexp_parser ():
    """
    Returns a  simple parser for  nested lists of real  numbers. Round
    parens () are assumed as customary in lisps.
    """

    # Punctuation literals (note round parens):
    LPAR, RPAR = map (Suppress, "()")

    # Real numbers:
    real_string = Regex (r"[+-]?\d+\.\d*([eE][+-]?\d+)?")
    real = real_string.setParseAction (lambda tokens: float (tokens[0]))

    # Voodoo:
    sexp = Forward ()
    sexp_list = Group (LPAR + ZeroOrMore (sexp) + RPAR)
    sexp << (real | sexp_list)

    return lambda s: sexp.parseString (s)[0]

from_sexp = make_sexp_parser ()


#
# Request  is  an s-expression  in  text  form  encoding an  array  of
# numbers, reply is an s-expression (e .  g) where g is itself a n x 3
# nested list  of numbers  again all in  text form.  The  parser turns
# nested list back into Python nested lists.
#
def protocol (server):
    """
    The result  is a  Func() derived from  the protocol  function fg()
    that implements  a "Taylor  expansion" of a  scalar PES at  x.  It
    returns a tuple of energy and gradients:

      fg: array (n x 3) -> number, array (n x 3)

    Server is supposed to turn text into text (meanigfully):

      server: text -> text

    For the properties  of Func() see module pts.func.   Ah, and there
    is another gotcha --- PTS python code operates in eV and angstroms
    whereas the server accepts angstroms but returns kcals.
    """

    def fg (x):
        # The parser  returns an object almost isomorphic  to a nested
        # list, but  not quite. In particular  array() constructor may
        # choke on such an object:
        y = from_sexp (server (to_sexp (x)))
        # print "XXX:\n", x, "->\n", y

        # Energy:
        e = y[0]

        # Gradients
        g = array (y[1:])

        return e * kcal, g * kcal

    #
    # The value  returned from here  is the value yielded  in Server()
    # context manager and is used in PTS code as this:
    #
    #   with Server (cmd) as f:
    #       # ... assuming f is a Func() here.
    #
    # So we have  to turn tuple valued fg()  into proper Func() either
    # here or  in the body of the  Server(). But since this  is also a
    # protocol thing we do it here:
    #
    return Func (taylor=fg)


#
# The bulk of python code  operates with callable Objects of type Func
# for historical  reasons. So that also a  server-based solution needs
# to  be cast  into that  pattern.  This  is supposed  to be  one such
# Func. Each  Func is also a context  manager and this is  a rare case
# when a Func implements a non-trivial context manager protocol.
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
    # subcommand   and   two  more   positional   arguments  will   be
    # appended. The  "server" subcommand must be  the first positional
    # argument in order to be interpreted as such:
    proc = Popen (args + ["server", inp, out])

    # This  function  takes  input  text  and returns  output  of  the
    # subprocess, both communicated via the fifos:
    def server (x):
        try:
            with open (inp, "w") as f:
                f.write (x)
            with open (out, "r") as f:
                y = f.read ()
        except KeyboardInterrupt as e:
            #
            # FIXME:  If the IO  is screwed,  we will  not be  able to
            # regularly shut  down the daemon.  Kill  it with SIGTERM.
            # Some minimal testing showed that proc.poll() will return
            # something  trueish  afterwards.   Also  if you  catch  a
            # generic  Exception here  the chances  are slim  that you
            # will catch  an asynchrounous KeyboardInterrupt  here and
            # now (see StackOverflow).
            #
            proc.terminate()
            proc.wait()
            raise e
        return y

    # This is a valid Func:
    func = protocol (server)

    # If an  exception happens in the  body of the  "with" statment we
    # will note it here:
    try:
        #
        # The control returns to the
        #
        #   with Server(...) as func:
        #
        # statement body:
        #
        yield func
    finally:
        #
        # When  leaving  the   "with"  statement  body  execute  this.
        # Executed  in  any  case,   even  if  exceptions  such  as  a
        # KeyboardInterrupt occur.  We do  not want zombie daemons and
        # stray files.
        #
        # Tell the daemon process to terminate.  This will not work if
        # user hit C-c while the  server() was already doing IO on the
        # fifos. In this case  the process should have been terminated
        # there.
        #
        if not proc.poll():
            # Regular shutdown sequence:
            with open (inp, "w") as f:
                f.write ("#f")  # convention
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
