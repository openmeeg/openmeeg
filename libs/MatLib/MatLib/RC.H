//  This code is from the book "More Effective C++" by Scott Meyers (1996).

#ifndef RC_H
#define RC_H

namespace utils {

    //	Template class for smart pointers-to-T objects;
    //	T must inherit from RCObject.

    template <class T>
    class RCPtr {
    public:
        
        RCPtr(T* realPtr = 0):  pointee(realPtr)     { init(); }
        RCPtr(const RCPtr& rhs):pointee(rhs.pointee) { init(); }

        ~RCPtr() { if (pointee) pointee->removeReference(); }

        RCPtr& operator=(const RCPtr&);

        T* operator->() const { return pointee;  }
        T& operator*()  const { return *pointee; }

    private:

        T* pointee;

        void init();
    };

    template <class T>
    void RCPtr<T>::init() {

        if (pointee==0) return;

        if (pointee->isShareable()==false)
            pointee = new T(*pointee);

        pointee->addReference();
    }

    template <class T>
    RCPtr<T>&
    RCPtr<T>::operator=(const RCPtr& rhs) {

        if (pointee!=rhs.pointee) {
            if (pointee) pointee->removeReference();

            pointee = rhs.pointee;
            init();
        }

        return *this;
    }

    //	Base class for reference-counted objects.

    class OPENMEEGMATHS_EXPORT RCObject {
    public:

        void addReference()      { ++refCount;                     }
        void removeReference()   { if (--refCount==0) delete this; }

        void markUnsharable()    { shareable = false;              }
        bool isShareable() const { return shareable;               }

        bool isShared() const    { return refCount>1;              }

    protected:

        RCObject(): refCount(0),shareable(true) { }

        RCObject(const RCObject&): refCount(0),shareable(true) { }

        virtual ~RCObject() { }
        
        RCObject& operator=(const RCObject&) { return *this; }

    private:

        mutable int    refCount;
        mutable bool   shareable;
    };

    template <class T>
    class RCIPtr {
    public:

        RCIPtr(T* realPtr=0):counter(new CountHolder) {
            counter->pointee = realPtr;
            init();
        }

        RCIPtr(const RCIPtr& rhs): counter(rhs.counter) { init(); }

        ~RCIPtr() { if (counter) counter->removeReference(); }

        RCIPtr& operator=(const RCIPtr&);

        T* operator->() const { return counter->pointee;    }
        T& operator*()  const { return *(counter->pointee); }

    private:

        struct CountHolder: public RCObject {
            ~CountHolder() { delete pointee; }

            T *pointee;
        };

        CountHolder *counter;

        void init();
    };

    template <class T>
    void RCIPtr<T>::init()
    {
        if (counter==0) return;

        if (counter->isShareable()==false) {
            counter = new CountHolder;
            counter->pointee = new T(*counter->pointee);
        }

        counter->addReference();
    }

    template <class T>
    RCIPtr<T>&
    RCIPtr<T>::operator=(const RCIPtr& rhs)
    {
        if (counter!=rhs.counter) {
            if (counter) counter->removeReference();
            counter = rhs.counter;
            init();
        }

        return *this;
    }
}

#endif  //  ! RC_H
