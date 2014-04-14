'use strict';

describe('Controller: TestTokenizerCtrl', function () {

  // load the controller's module
  beforeEach(module('clothoApp'));

  var TestTokenizerCtrl,
    scope;

  // Initialize the controller and a mock scope
  beforeEach(inject(function ($controller, $rootScope) {
    scope = $rootScope.$new();
    TestTokenizerCtrl = $controller('TestTokenizerCtrl', {
      $scope: scope
    });
  }));

  it('should attach a list of awesomeThings to the scope', function () {
    expect(scope.awesomeThings.length).toBe(3);
  });
});
