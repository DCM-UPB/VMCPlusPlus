namespace std {
  // STL vector class
  template<class T> class vector { public: T element; };

  // STL array class
  template<class T, size_t N> class array { public: T element; public: N size; };

  // STL set class
  template<class T> class set { public: T element; };

  // STL map class
  template<class K, class V> class map { public: K Key; public: V Value; };

  // STL pair class
  template<class A, class B> class pair { public: A first; public: B first; };

  // unique_ptr
  template<class T> class unique_ptr { public: *T rawptr; };
}
