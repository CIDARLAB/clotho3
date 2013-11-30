'use strict';

describe('Directive: sharable', function () {

  // load the directive's module
  beforeEach(module('clothoApp'));

  var element,
    scope;

  beforeEach(inject(function ($rootScope) {
    scope = $rootScope.$new();
  }));

  it('should make hidden element visible', inject(function ($compile) {
    element = angular.element('<sharable></sharable>');
    element = $compile(element)(scope);
    expect(element.text()).toBe('this is the sharable directive');
  }));
});
