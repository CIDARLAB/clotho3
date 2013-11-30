'use strict';

describe('Directive: darkBackground', function () {

  // load the directive's module
  beforeEach(module('clothoApp'));

  var element,
    scope;

  beforeEach(inject(function ($rootScope) {
    scope = $rootScope.$new();
  }));

  it('should make hidden element visible', inject(function ($compile) {
    element = angular.element('<dark-background></dark-background>');
    element = $compile(element)(scope);
    expect(element.text()).toBe('this is the darkBackground directive');
  }));
});
