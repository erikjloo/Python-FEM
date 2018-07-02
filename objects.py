class Object(object):
    """ Object class """

    def accepts(*types):
        def check_accepts(f):
            assert len(types) == f.func_code.co_argcount

            def new_f(*args, **kwds):
                for (a, t) in zip(args, types):
                    assert isinstance(a, t), \
                        "arg %r does not match %s" % (a, t)
                return f(*args, **kwds)
            new_f.func_name = f.func_name
            return new_f
        return check_accepts

    def enforce(*types):
        def decorator(f):
            def new_f(*args, **kwds):
                #we need to convert args into something mutable
                newargs = []
                for (a, t) in zip(args, types):
                    # feel free to have more elaborated convertion
                    newargs.append(t(a))
                return f(*newargs, **kwds)
            return new_f
        return decorator