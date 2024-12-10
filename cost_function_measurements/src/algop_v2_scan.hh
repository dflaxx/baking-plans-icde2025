#pragma once

#include "algop_v2_base.hh"
#include "RelationRS.hh"


/*
 * Table Scan Operator
 */
template <alg2_consumer_c Tconsumer> 
class AlgOpScan2 : public AlgOpBase2 {
  public:
    using consumer_t  = Tconsumer;
    using input_t     = typename consumer_t::input_t;  // push input as is
    using output_t    = typename consumer_t::input_t;
    using input_rel_t = RelationRS2<input_t>;
  public:
    inline AlgOpScan2(consumer_t* aConsumer)
      : AlgOpBase2(algop_et::kAlgOpScan), _consumer(aConsumer) {}
  public:
    void run(const input_rel_t& aRelation);
  public:
    inline bool mem_alloc() { assert(nullptr != _consumer); return _consumer->mem_alloc(); }
    inline bool mem_init() { assert(nullptr != _consumer); return _consumer->mem_init(); }
    inline void clear() { assert(nullptr != _consumer); _consumer->clear(); }
    inline void mem_free() { assert(nullptr != _consumer); _consumer->mem_free(); }
  public:
    inline const consumer_t* consumer() const { return _consumer; }
    inline const AlgOpBase2* consumer_poly() const { return _consumer; }
  private:
    consumer_t*  _consumer;
};

template <alg2_consumer_c Tconsumer>
void
AlgOpScan2<Tconsumer>::run(const input_rel_t& aRelation) {
  reset();
  _consumer->init();

  assert(0 == count());
  if (nullptr != consumer()) {
    assert(0 == consumer()->count());
    if (nullptr != consumer()->consumer()) {
      assert(0 == consumer()->consumer_poly()->count());
    }
  }

  for (size_t i = 0; i < aRelation.card(); ++i) {
    inc();
    _consumer->step(aRelation.ptr(i));
  }
  _consumer->fin();
}
