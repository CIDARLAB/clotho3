'use strict';

describe('Service: EditorSchemas', function () {

  // load the service's module
  beforeEach(module('clothoApp'));

  // instantiate service
  var EditorSchemas;
  beforeEach(inject(function (_EditorSchemas_) {
    EditorSchemas = _EditorSchemas_;
  }));

  it('should do something', function () {
    expect(!!EditorSchemas).toBe(true);
  });

});
